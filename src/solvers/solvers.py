#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue Jun  7 13:18:01 2022

@author: sampasmann
"""
import sys
sys.path.append("/Users/sampasmann/Documents/GitHub/QMC1D/")
import numpy as np
from src.solvers.maps import SI_Map, RHS, MatVec_data, MatVec
from scipy.sparse.linalg import gmres, bicgstab, LinearOperator



def Picard(qmc_data,tol=1e-5,maxit=50,report_progress=False):
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
    if (report_progress):
        print("--------- Source Iteration ---------")
        print("Material: ", qmc_data.material_code)
        print("Random Number Generator: ", qmc_data.generator)
        print("Number of Particles per Iteration: ", qmc_data.N)
        print("Number of Spatial Cells: ", qmc_data.Nx)
    Nx      = qmc_data.Nx
    G       = qmc_data.G
    phi0    = qmc_data.source
    itc     = 0
    diff    = 1.0
    phic    = np.copy(phi0)
    phi     = np.copy(phi0)
    reshist = np.empty(0)
    
    while (itc < maxit) and (diff > tol):
        phi_out = SI_Map(phic, qmc_data)
        phi     = np.reshape(phi_out,(Nx,G))
        diff    = np.linalg.norm((phic-phi))
        reshist = np.append(reshist, diff)
        phic    = np.copy(phi)
        itc += 1
        if (report_progress):
            print("**********************")
            print("Iteration:", itc, "change: ",diff)
    #tallies = out[1]
    return phi



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
    Nx       = qmc_data.Nx
    G        = qmc_data.G
    Nv       = Nx*G
    global itt
    itt = 0
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

    gmres_out = gmres(A,b,x0=phi0,tol=tol,maxiter=maxit)
    phi = gmres_out[0]
    phi = np.reshape(phi, (Nx,G))
    
    if (gmres_out[1]>0):
        print("Convergence to tolerance not achieved: Maximum number of iterations.")
    elif (gmres_out[1]<0):
        print("Illegal input or breakdown")
        
    return phi

def BICGSTAB(qmc_data,tol=1e-5,maxit=50):
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
    Nx       = qmc_data.Nx
    G        = qmc_data.G
    Nv       = Nx*G
    global itt
    itt = 0
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
    
    if (gmres_out[1]>0):
        print("Convergence to tolerance not achieved: Maximum number of iterations.")
    elif (gmres_out[1]<0):
        print("Illegal input or breakdown")
        
    return phi


if (__name__ == "__main__"):
    from src.init_files.mg_init import MultiGroupInit, TrueFlux
    import time 
    
    N = 2**10
    Nx = 10
    maxit = 20
    tol = 1e-5
    qmc_data = MultiGroupInit(N=N, Nx=Nx, generator="halton")
    
    start = time.time()
    phi_gmres = GMRES(qmc_data,tol=tol,maxit=maxit)
    stop = time.time()
    print("GMRES took: ", stop-start, " seconds")
    
    start = time.time()
    phi_bic = BICGSTAB(qmc_data,tol=tol,maxit=maxit)
    stop = time.time()
    print("BiCGSTAB took: ", stop-start, " seconds")
    
    start = time.time()
    phi_pic = Picard(qmc_data,tol=tol,maxit=maxit)
    stop = time.time()
    print("Picard took: ", stop-start, " seconds")
    
    Phi_sol = TrueFlux(qmc_data.material, qmc_data.source, Nx)
