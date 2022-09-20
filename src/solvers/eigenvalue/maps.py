#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue Jun  7 13:17:55 2022

@author: sampasmann
"""
import numpy as np
from src.functions.tallies import Tallies
from src.functions.sweep import Sweep
from src.functions.source import GetCriticalitySource
from mpi4py import MPI


def SI_Map(phi_f, phi_in, qmc_data):
    """
    PI_Map(phi_in, qmc_data)
    -----------------------
    Source Iteration Map
    """
    G   = qmc_data.G
    Nx  = qmc_data.Nx
    Nv  = phi_in.size
    try:
        (Nv == Nx*G)  
    except Exception as e: print(e) 
    phi_in      = np.reshape(phi_in, (Nx,G))
    
    tallies     = Tallies(qmc_data)
    source      = GetCriticalitySource(phi_f, phi_in, qmc_data)
    sweep       = Sweep(qmc_data) # samples are gneratred with initialization of sweep
    #print("         QMC Sweep ")
    sweep.Run(tallies, source) # QMC sweep
    phi_out     = tallies.phi_avg
    phi_out     = np.reshape(phi_out,(Nv,1))
    # all reduce phi_out here (they automatically wait for each other)
    comm        = MPI.COMM_WORLD
    phi_out     = comm.allreduce(phi_out,op=MPI.SUM)
    
    return phi_out


def RHS(qmc_data):
    """
    RHS(qmc_data)
    -------------
    We solve A x = b with a Krylov method. This function extracts
    b from Sam's qmc_data structure by doing a transport sweep with
    zero scattering term.
    """
    G       = qmc_data.G
    Nx      = qmc_data.Nx
    Nv      = Nx*G
    phi_f   = qmc_data.phi_f
    zed     = np.zeros((Nx,G))
    b       = SI_Map(phi_f, zed, qmc_data) # qmc_sweep with phi(0)
    
    return b


def MatVec_data(qmc_data):
    """
    MXV_data(qmc_data)
    ------------------
    This function adds the right side of the linear system to Sam's 
    qmc_data structure so I can pass it to the matrix-vector product.
    """
    b           = RHS(qmc_data)
    global matvec_data
    matvec_data = [b, qmc_data]
    return matvec_data


def MatVec(phi_in):
    """
    MXV(phi_in)
    ---------------------
    We solve A x = b with a Krylov method. This function extracts the
    matrix-vector product A * phi_in from Sam's qmc_data structure by 
    doing a transport sweep with zero boundary conditions and zero external
    source.
    """
    b           = matvec_data[0]
    qmc_data    = matvec_data[1]
    phi_f       = qmc_data.phi_f
    Nx          = qmc_data.Nx
    G           = qmc_data.G
    Nv          = Nx*G
    phi_in      = np.reshape(phi_in,(Nv,1))

    #qmc_data.source = np.zeros((Nx,G))
    axv             = phi_in - SI_Map(phi_f, phi_in, qmc_data) + b

    return axv
