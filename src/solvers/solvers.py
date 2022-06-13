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
from scipy.sparse.linalg import gmres, LinearOperator

from src.init_files.mg_init import MultiGroupInit
from src.init_files.garcia_init import GarciaInit
#from src.functions.save_data import SaveData


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
    return [phi, reshist]



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
    matvec_data = MatVec_data(qmc_data)
    A        = LinearOperator((Nx,G), matvec=MatVec) # this line is the problem
    b        = matvec_data[0]
    phi0     = qmc_data.source
    
    gmres_out = gmres(A,b,x0=phi0,tol=tol,maxiter=maxit)
    phi = gmres_out[0]
    
    return phi



if (__name__ == "__main__"):
    Nx = 10
    qmc_data = MultiGroupInit(Nx=Nx, generator="halton")
    maxit = 20
    phi_out = GMRES(qmc_data,maxit=maxit)
    #out = Picard(qmc_data,maxit=maxit,report_progress=True)
    #tallies = out[2]