#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue Jun  7 13:18:01 2022

@author: sampasmann
"""
import numpy as np
import time
from src.functions.tallies import Tallies
from src.functions.sweep import Sweep
from src.functions.source import GetCriticalitySource
from src.functions.save_data import SaveData
from src.solvers.eigenvalue.maps import SI_Map, RHS, MatVec_data, MatVec
from scipy.sparse.linalg import gmres, lgmres, bicgstab, LinearOperator
from mpi4py import MPI


def PowerIteration(qmc_data):
    solver          = "LGMRES"
    itt             = 0
    max_SI_iter     = 10
    max_PI_iter     = 10
    k_tol           = 1e-3
    phi_tol         = 1e-1
    k               = qmc_data.keff
    dk              = 1.0
    phi_old         = qmc_data.phi_f.copy()
    
    print("--------- Power Iteration ---------")
    print("Inner Sovler:",                       solver)
    print("Material: ",                          qmc_data.material_code)
    print("Random Number Generator: ",           qmc_data.generator)
    print("Number of Particles per Iteration: ", qmc_data.N)
    print("Number of Spatial Cells: ",           qmc_data.Nx)
    print("Initial K: ",                         qmc_data.keff)
    
    # iterate over k effective
    while (itt<=max_PI_iter) and (dk>=k_tol):
        # iterate over scattering source
        phi_new         = InnerIteration(qmc_data, solver=solver, maxit=max_SI_iter,tol=phi_tol)
        k_old           = k
        k               = UpdateK(phi_old, phi_new, qmc_data)
        qmc_data.keff   = k
        qmc_data.phi_f  = phi_new.copy()
        phi_old         = phi_new.copy()
        dk              = abs(k-k_old)
        itt             += 1
        
        print("**********************")
        print("Iteration:", itt, "dk: ",dk)
        print("k: ", k)
    if (itt>=max_PI_iter):
        print("Power Iteration convergence to tolerance not achieved: Maximum number of iterations.")
    elif (dk<=k_tol):
        print("-------------------------------")
        print("Successful Power Iteration convergence.")

def UpdateK(phi_f, phi_s, qmc_data):
    keff        = qmc_data.keff
    material    = qmc_data.material
    keff        *= (np.sum(material.nu*material.sigf*phi_s)
                    /np.sum(material.nu*material.sigf*phi_f))
    return keff


def InnerIteration(qmc_data,solver="LGMRES",tol=1e-5,maxit=50,save_data=False):
    """
    Parameters
    ----------
    qmc_data : TYPE
        DESCRIPTION.
    tol : TYPE, optional
        DESCRIPTION. The default is 1e-5.
    maxit : TYPE, optional
        DESCRIPTION. The default is 50.

    Returns
    -------
    phi : TYPE
        DESCRIPTION.

    """
    comm        = MPI.COMM_WORLD
    rank        = comm.Get_rank()
    nproc       = comm.Get_size()
    Nx          = qmc_data.Nx
    G           = qmc_data.G
    Nv          = Nx*G
    start       = time.time()
    matvec_data = MatVec_data(qmc_data)
    A           = LinearOperator((Nv,Nv), 
                              matvec=MatVec,
                              rmatvec=MatVec,
                              matmat= MatVec,
                              rmatmat=MatVec,
                              dtype=float)
    b           = matvec_data[0]
    phi0        = qmc_data.source
    phi0        = np.reshape(phi0,(Nv,1))

    print("     Inner Iteration: ")
    if (solver=="LGMRES"):
        gmres_out   = lgmres(A,b,x0=phi0,tol=tol,maxiter=maxit)
    elif (solver=="GMRES"):
        gmres_out   = gmres(A,b,x0=phi0,tol=tol,maxiter=maxit)
    elif (solver=="BICGSTAB"):
        gmres_out   = bicgstab(A,b,x0=phi0,tol=tol,maxiter=maxit)
    else:
        print(" Not a valid solver ")
        Exception
        
    phi         = gmres_out[0]
    exitCode    = gmres_out[1]
    phi         = np.reshape(phi, (Nx,G))
    stop        = time.time()
    run_time    = stop - start
    
    if (rank==0):
        if (save_data):
            sim_data = SimData(phi, run_time, tol, nproc)
            SaveData(qmc_data, sim_data)
        if (exitCode>0):
            print("     Convergence to tolerance not achieved: Maximum number of iterations.")
        elif (exitCode<0):
            print("     Illegal input or breakdown.")
        elif (exitCode==0):
            print("     Successful convergence.")
        
    return phi


def SimData(phi, time, tol, nproc):
    data = {
        "phi": phi,
        "run_time": time,
        "tolerance": tol,
        "nproc": nproc
        }
    return data


