#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue Jun  7 13:17:55 2022

@author: sampasmann
"""
import numpy as np
from src.functions.tallies import Tallies
from src.functions.sweep import Sweep
from src.functions.source import GetSource

def SI_Map(phi_in, qmc_data):
    """
    SI_Map(phi_in, qmc_data)
    -----------------------
    Source Iteration Map
    """
    G   = qmc_data.G
    Nx  = qmc_data.Nx
    Nv  = phi_in.size
    try:
        (Nv == Nx*G)  
    except Exception as e: print(e) 
    phi_avg = np.reshape(phi_in, (Nx,G))
    source  = GetSource(phi_avg, qmc_data)
    tallies = Tallies(qmc_data)
    sweep   = Sweep(qmc_data)
    sweep.Run(tallies, source)
    phi_out = tallies.phi_avg
    phi_out = np.reshape(phi_out,Nv,)
    return phi_out


def RHS(qmc_data):
    """
    RHS(qmc_data)
    -------------
    We solve A x = b with a Krylov method. This function extracts
    b from Sam's qmc_data structure by doing a transport sweep with
    zero scattering term.
    """
    G   = qmc_data.G
    Nx  = qmc_data.Nx
    Nv  = Nx*G
    zed = np.zeros(Nv)
    bout = SI_Map(zed,qmc_data)
    return bout


def MXV_data():
    """
    MXV_data(qmc_data)
    ------------------
    This function adds the right side of the linear system to Sam's 
    qmc_data structure so I can pass it to the matrix-vector product.
    """
    b = RHS(qmc_data)
    mxv_data = (b, qmc_data)
    return mxv_data


def MXV(phi_in, mxv_data):
    """
    MXV(phi_in, mxv_data)
    ---------------------
    We solve A x = b with a Krylov method. This function extracts the
    matrix-vector product A * phi_in from Sam's qmc_data structure by 
    doing a transport sweep with zero boundary conditions and zero external
    source.
    """
    b           = mxv_data[0]
    qmc_data    = mxv_data[1]
    mxvp        = SI_Map(phi_in ,qmc_data)
    mxv         = mxvp - b
    axv         = phi_in - mxv
    return axv







