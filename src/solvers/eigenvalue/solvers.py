#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue Jun  7 13:18:01 2022

@author: sampasmann
"""

import numpy as np
import time
from mpi4py import MPI
from src.solvers.eigenvalue import PI_Map, RHS, MatVec_data, MatVec
from src.functions.save_data import SaveData
from scipy.sparse.linalg import gmres, lgmres, bicgstab, LinearOperator



def GMRES(qmc_data,tol=1e-5,maxit=50):
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
    comm = MPI.COMM_WORLD
    rank = comm.Get_rank()
    nproc = comm.Get_size()
    
    Nx        = qmc_data.Nx
    G         = qmc_data.G
    Nv        = Nx*G
    start = time.time()
    matvec_data = MatVec_data(qmc_data)
    A         = LinearOperator((Nv,Nv), 
                              matvec=MatVec,
                              rmatvec=MatVec,
                              matmat= MatVec,
                              rmatmat=MatVec,
                              dtype=float) # this line is the problem
    b         = matvec_data[0]
    phi0      = qmc_data.source
    phi0      = np.reshape(phi0,(Nv,1))

    gmres_out = gmres(A,b,x0=phi0,tol=tol,maxiter=maxit)
    phi       = gmres_out[0]
    phi       = np.reshape(phi, (Nx,G))
    stop      = time.time()
    run_time  = stop-start
    
    if (rank==0):
        sim_data = SimData(phi, run_time, tol, nproc)
        SaveData(qmc_data, sim_data)
        if (gmres_out[1]>0):
            print("Convergence to tolerance not achieved: Maximum number of iterations.")
        elif (gmres_out[1]<0):
            print("Illegal input or breakdown")
        
    return phi

def LGMRES(qmc_data,tol=1e-5,maxit=50,save_data=True):
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
    comm = MPI.COMM_WORLD
    rank = comm.Get_rank()
    nproc = comm.Get_size()
    
    Nx       = qmc_data.Nx
    G        = qmc_data.G
    Nv       = Nx*G
    start = time.time()
    matvec_data = MatVec_data(qmc_data)
    A        = LinearOperator((Nv,Nv), 
                              matvec=MatVec,
                              rmatvec=MatVec,
                              matmat= MatVec,
                              rmatmat=MatVec,
                              dtype=float)
    b        = matvec_data[0]
    phi0     = qmc_data.source
    phi0 = np.reshape(phi0,(Nv,1))

    gmres_out = lgmres(A,b,x0=phi0,atol=tol,maxiter=maxit)
    phi = gmres_out[0]
    phi = np.reshape(phi, (Nx,G))
    stop = time.time()
    run_time = stop - start
    
    if (rank==0):
        if (save_data):
            sim_data = SimData(phi, run_time, tol, nproc)
            SaveData(qmc_data, sim_data)
        if (gmres_out[1]>0):
            print("Convergence to tolerance not achieved: Maximum number of iterations.")
        elif (gmres_out[1]<0):
            print("Illegal input or breakdown")
        
    return phi

def BICGSTAB(qmc_data,tol=1e-5,maxit=50,save_data=True):
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
    comm = MPI.COMM_WORLD
    rank = comm.Get_rank()
    nproc = comm.Get_size()
    
    Nx       = qmc_data.Nx
    G        = qmc_data.G
    Nv       = Nx*G
    start = time.time()
    
    matvec_data = MatVec_data(qmc_data)
    A        = LinearOperator((Nv,Nv), 
                              matvec=MatVec,
                              rmatvec=MatVec,
                              matmat= MatVec,
                              rmatmat=MatVec,
                              dtype=float) # this line is the problem
    b        = matvec_data[0]
    phi0     = qmc_data.source
    phi0 = np.reshape(phi0,(Nv,1))

    gmres_out = bicgstab(A,b,x0=phi0,tol=tol,maxiter=maxit)
    phi = gmres_out[0]
    phi = np.reshape(phi, (Nx,G))
    stop = time.time()
    run_time = stop - start
    
    if (rank==0):
        if (save_data):
            sim_data = SimData(phi, run_time, tol, nproc)
            SaveData(qmc_data, sim_data)
        if (gmres_out[1]>0):
            print("Convergence to tolerance not achieved: Maximum number of iterations.")
        elif (gmres_out[1]<0):
            print("Illegal input or breakdown")
        
    return phi


def SimData(phi, time, tol, nproc):
    data = {
        "phi": phi,
        "run_time": time,
        "tolerance": tol,
        "nproc": nproc
        }
    return data


