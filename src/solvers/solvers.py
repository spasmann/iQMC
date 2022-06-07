#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue Jun  7 13:18:01 2022

@author: sampasmann
"""
import sys
sys.path.append("/Users/sampasmann/Documents/GitHub/QMC1D/")
import numpy as np
from src.solvers.maps import SI_Map, RHS, MXV_data, MXV
from scipy.sparse.linalg import gmres, LinearOperator

from src.init_files.mg_init import MultiGroupInit
#from src.functions.save_data import SaveData


def Picard(phi0,SI_Map,tol,maxit,qmc_data,report_progress=False):
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
    itc     = 0
    diff    = 1.0
    phic    = np.copy(phi0)
    phi     = np.copy(phi0)
    reshist = np.empty(0)
    
    while (itc < maxit) and (diff > tol):
        out     = SI_Map(phic, qmc_data)
        phi     = np.reshape(out[0],(Nx,G))
        diff    = np.linalg.norm((phic-phi))
        reshist = np.append(reshist, diff)
        phic    = np.copy(phi)
        itc += 1
        if (report_progress):
            print("**********************")
            print("Iteration:", itc, "change: ",diff)
    tallies = out[1]
    return [phi, reshist, tallies]



def GMRES(phi0,SI_Map,tol,maxit,qmc_data):
    """
    Parameters
    ----------
    phi0 : TYPE
        DESCRIPTION.
    SI_Map : TYPE
        DESCRIPTION.
    tol : TYPE
        DESCRIPTION.
    maxit : TYPE
        DESCRIPTION.
    qmc_data : TYPE
        DESCRIPTION.

    Returns
    -------
    None.
    """
    Nx       = qmc_data.Nx
    G        = qmc_data.G
    mxv_data = MXV_data(qmc_data)
    A        = LinearOperator((Nx,G), matvec=MXV)
    b        = mxv_data[0]
    phi0     = qmc_data.source
    
    gmres_out = gmres(A,b,x0=phi0,tol=tol,maxiter=maxit)
    phi = gmres_out[0]
    
    return phi



if (__name__ == "__main__"):
    
    qmc_data = MultiGroupInit(generator="halton")
    tol = 1e-5
    maxit = 25
    phi0 = qmc_data.source
    phi_out = GMRES(phi0,SI_Map,tol,maxit,qmc_data)
    #out = Picard(phi0,SI_Map,tol,maxit,qmc_data,report_progress=True)
    #tallies = out[2]