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

from src.init_files.garcia_init import GarciaInit


def Picard(phi0,SI_Map,tol,maxit,qmc_data,report_progress=False):
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
    itc     = 0
    diff    = 1.0
    phic    = np.copy(phi0)
    phi     = np.copy(phi0)
    reshist = np.empty(0)
    
    while (itc < maxit) and (diff > tol):
        phi     = SI_Map(phic, qmc_data)
        diff    = np.linalg.norm((phic-phi))
        reshist = np.append(reshist, diff)
        phic    = np.copy(phi)
        itc += 1
        if (report_progress):
            print("**********************")
            print("Iteration:", itc, "change: ",diff)
    return [phi, reshist]



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
    
    return [phi, reshist]



if (__name__ == "__main__"):
    qmc_data = GarciaInit(generator="halton")
    tol = 1e-5
    maxit = 25
    phi0 = qmc_data.source
    Picard(phi0,SI_Map,tol,maxit,qmc_data,report_progress=True)