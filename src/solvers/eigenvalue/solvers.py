#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue Jun  7 13:18:01 2022

@author: sampasmann
"""
import time
import numpy as np
from mpi4py import MPI
from src.functions.save_data import SaveData
from src.solvers.fixed_source.solvers import Picard
from src.solvers.eigenvalue.maps import MatVec_data, MatVec
from scipy.sparse.linalg import gmres, lgmres, bicgstab, LinearOperator


class gmres_counter(object):
    def __init__(self, disp=True):
        self._disp = disp
        self.iter = 0
        self.callbacks = []
    def __call__(self, rk=None):
        self.callbacks.append(rk.copy())
        self.iter += 1
        if self._disp:
            #print(rk)
            if (self.iter>1):
                print("     Iteration:", self.iter-1, "change: ",np.linalg.norm((rk - self.callbacks[self.iter-2])))
                
def PowerIteration(qmc_data, solver="LGMRES", max_outter_itt=10, max_inner_itt=10, outter_tol=1e-5, inner_tol=1e-5):
    comm        = MPI.COMM_WORLD
    rank        = comm.Get_rank()
    nproc       = comm.Get_size()
    
    itt             = 0
    k               = qmc_data.keff
    dk              = 1.0
    phi_old         = qmc_data.phi_f.copy()
    #res_hist        = []
    k_hist          = []
    if (rank==0):
        print("--------- K-Effective Eigenvalue Problem ---------")
        print("Outter Solver: Power Iteration")
        print("Inner Sovler:",                       solver)
        print("Material: ",                          qmc_data.material_code)
        print("Random Number Generator: ",           qmc_data.generator)
        print("Number of Particles per Iteration: ", qmc_data.N)
        print("Number of Spatial Cells: ",           qmc_data.Nx)
        print("Initial K: ",                         qmc_data.keff)
    
    # iterate over k effective
    while (itt<=max_outter_itt) and (dk>=outter_tol):
        # iterate over scattering source
        phi_new         = InnerIteration(qmc_data, solver=solver, maxit=max_inner_itt,tol=inner_tol)
        #phi_hist.append(phi_new)
        k_old           = k
        k               = UpdateK(phi_old, phi_new, qmc_data)
        k_hist.append(k)
        qmc_data.keff   = k
        #res_hist.append(np.linalg.norm(phi_new-phi_old))
        qmc_data.phi_f  = phi_new.copy()
        phi_old         = phi_new.copy() # /norm(phi_new)
        dk              = abs(k-k_old)
        itt             += 1
        if (rank==0):
            print("**********************")
            print("Iteration:", itt)
            print("k: ", k)
            print("dk: ",dk)
    if (rank==0):
        if (itt>=max_outter_itt):
            print("Power Iteration convergence to tolerance not achieved: Maximum number of iterations.")
        elif (dk<=outter_tol):
            print("-------------------------------")
            print("Successful Power Iteration convergence.")
            
    return phi_new, k_hist #, res_hist

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

    if (rank==0):
        print("     Inner Iteration: ")
    if (solver=="LGMRES"):
        counter     = gmres_counter()
        gmres_out   = lgmres(A,b,x0=phi0,tol=tol,maxiter=maxit, callback=counter)
    elif (solver=="GMRES"):
        counter     = gmres_counter()
        gmres_out   = gmres(A,b,x0=phi0,tol=tol,maxiter=maxit, callback=counter)
    elif (solver=="BICGSTAB"):
        counter     = gmres_counter()
        gmres_out   = bicgstab(A,b,x0=phi0,tol=tol,maxiter=maxit, callback=counter)
    elif (solver=="Picard"):
        gmres_out = Picard(qmc_data,tol=tol,maxit=maxit,save_data=False,report_progress=True)
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


