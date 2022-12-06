#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue Jun  7 13:18:01 2022

@author: sampasmann
"""

import numpy as np
import time
from mpi4py import MPI
from src.solvers.fixed_source.maps import SI_Map, RHS, MatVec_data, MatVec
from src.functions.save_data import SaveData
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
                print("**********************")
                print("Iteration:", self.iter-1, "change: ",np.linalg.norm((rk - self.callbacks[self.iter-2])))
            
            
def Picard(qmc_data,tol=1e-5,maxit=40,save_data=False,report_progress=True):
    """
    Parameters
    ----------
    phi0 : Starting volumetric source
    SI_Map : Source Iteration Map Function
    tol : Iteration tolerance
    maxit : Maximum number of iterations
    qmc_data : Object from init_files
    report_progress: boolean, print progress of iterative method

    Returns
    -------
    None.

    """
    comm    = MPI.COMM_WORLD
    rank    = comm.Get_rank()
    nproc   = comm.Get_size()

    if (qmc_data.source_tilt):
        phi0 = np.append(qmc_data.tallies.phi_avg, qmc_data.tallies.dphi_s)
    else:
        phi0 = qmc_data.tallies.phi_avg
    itc     = 0
    diff    = 1.0
    phic    = np.copy(phi0)
    reshist = np.empty(0)
    start   = time.time()
    while (itc < maxit) and (diff > tol):
        phi_out = SI_Map(phic, qmc_data)
        diff    = np.linalg.norm((phic-phi_out))
        reshist = np.append(reshist, diff)
        phic    = np.copy(phi_out)
        itc += 1
        if (report_progress) and (rank==0):
            print("**********************")
            print("Iteration:", itc, "change: ",diff)
    stop = time.time()
    run_time = stop-start
    if (rank==0):
        if (save_data==True):
            sim_data = SimData(phi_out, run_time, tol, nproc)
            SaveData(qmc_data, sim_data)
            
    if (qmc_data.source_tilt):
        Nv = int(qmc_data.Nx*qmc_data.G)
        phi_out = phi_out[:Nv]
    
    return phi_out


def FixedSource(qmc_data, solver="LGMRES", tol=1e-5, maxit=100, report_progress=True, save_data=False):
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
    start       = time.time()
    Nx          = qmc_data.Nx
    G           = qmc_data.G
    
    if (report_progress) and (rank==0):
        print("--------- Fixed Source Problem ---------")
        print("Solver:",                             solver)
        print("Material: ",                          qmc_data.material_code)
        print("Random Number Generator: ",           qmc_data.generator)
        print("Number of Particles per Iteration: ", qmc_data.N)
        print("Number of Spatial Cells: ",           qmc_data.Nx)

    if (solver == "Picard"):
        phi         = Picard(qmc_data,maxit=maxit,tol=tol,report_progress=report_progress)
        exitCode    = 0
    else:
        Nt          = qmc_data.Nt
        matvec_data = MatVec_data(qmc_data)
        A           = LinearOperator((Nt,Nt), 
                                  matvec=MatVec,
                                  rmatvec=MatVec,
                                  matmat= MatVec,
                                  rmatmat=MatVec,
                                  dtype=float)
        b           = matvec_data[0]
        
        if (qmc_data.source_tilt):
            phi0 = np.append(qmc_data.tallies.phi_avg, qmc_data.tallies.dphi_s)
        else:
            phi0 = qmc_data.tallies.phi_avg
        phi0 = np.reshape(phi0,(Nt,1))
        if (solver=="LGMRES"):
            counter     = gmres_counter(disp=report_progress)
            gmres_out   = lgmres(A,b,x0=phi0,tol=tol,maxiter=maxit,callback=counter)
            phi         = gmres_out[0]
            exitCode    = gmres_out[1]
        elif (solver=="GMRES"):
            counter     = gmres_counter(disp=report_progress)
            gmres_out   = gmres(A,b,x0=phi0,tol=tol,maxiter=maxit,callback=counter)
            phi         = gmres_out[0]
            exitCode    = gmres_out[1]
        elif (solver=="BICGSTAB"):
            counter     = gmres_counter(disp=report_progress)
            gmres_out   = bicgstab(A,b,x0=phi0,tol=tol,maxiter=maxit,callback=counter)
            phi         = gmres_out[0]
            exitCode    = gmres_out[1]
        else:
            if (rank==0):
                print(" Invalid Solver ")
                Exception
                
    phi         = np.reshape(phi[:int(Nx*G)], (Nx,G))        
    stop        = time.time()
    run_time    = stop - start
    
    if (rank==0):
        if (save_data):
            sim_data = SimData(phi, run_time, tol, nproc)
            SaveData(qmc_data, sim_data)
        if (exitCode>0) and (report_progress):
            print("Convergence to tolerance not achieved: Maximum number of iterations.")
        elif (exitCode<0) and (report_progress):
            print("Illegal input or breakdown")
        elif (exitCode==0) and (report_progress):
            print("Successful Convergence.")
        
    return phi


def SimData(phi, time, tol, nproc):
    data = {
        "phi": phi,
        "run_time": time,
        "tolerance": tol,
        "nproc": nproc
        }
    return data

    
