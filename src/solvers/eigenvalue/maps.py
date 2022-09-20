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
from src.solvers.eigenvalue.solvers import LGMRES, Inner_LGMRES
from mpi4py import MPI


def UpdateK(phi_f, phi_s, qmc_data):
    keff        = qmc_data.keff
    material    = qmc_data.material
    keff        *= (np.sum(material.nu*material.sigf*phi_s)
                    /np.sum(material.nu*material.sigf*phi_f))
    return keff
        
        
def SI_Map(phi_f, phi_s, qmc_data):
    """
    Parameters
    ----------
    phi_in_s : scalar flux for inner scatter source iteration
    phi_in_f : scalar flux for outter fission source iteration
    qmc_data : TYPE DESCRIPTION.

    Returns
    -------
    phi_out : 

    """

    G   = qmc_data.G
    Nx  = qmc_data.Nx
    Nv  = phi_in.size
    try:
        (Nv == Nx*G)  
    except Exception as e: print(e) 
    phi_s = np.reshape(phi_in_s, (Nx,G))
    phi_f = np.reshape(phi_in_f, (Nx,G))
    
    source  = GetCriticalitySource(keff, phi_s, phi_f, qmc_data)
    tallies = Tallies(qmc_data)
    sweep   = Sweep(qmc_data) # samples are gneratred with initialization of sweep
    sweep.Run(tallies, source) # QMC sweep
    phi_out = tallies.phi_avg
    phi_out = np.reshape(phi_out,(Nv,1))
    
    # all reduce phi_out here (processors automatically wait for each other)
    comm    = MPI.COMM_WORLD
    phi_out = comm.allreduce(phi_out,op=MPI.SUM)
    
    return phi_out


def PI_Map(phi_in, qmc_data):
    # phi_f == phi_in
    # phi_s == phi_out
    G               = qmc_data.G
    Nx              = qmc_data.Nx
    Nv              = phi_in.size
    qmc_data.source = phi_in
    phi_out         = Inner_LGMRES(qmc_data,tol=1e-5,maxit=50,save_data=False)
    phi_out         = np.reshape(phi_out,(Nv,1))
    
    qmc_data.k      = UpdateK(phi_in, phi_out, qmc_data)
    qmc_data.phi_f  = phi_out
    
    return phi_out


def PI_RHS(qmc_data):
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
    zed = np.zeros((Nx,G))
    b   = SI_Map(zed, zed, qmc_data) # qmc_sweep with phi(0)
    return b

def SI_RHS(qmc_data):
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
    zed = np.zeros((Nx,G))
    b   = SI_Map(zed, zed, qmc_data) # qmc_sweep with phi(0)
    return b


def PI_MatVec_data(qmc_data):
    """
    MXV_data(qmc_data)
    ------------------
    This function adds the right side of the linear system to Sam's 
    qmc_data structure so I can pass it to the matrix-vector product.
    """
    b   = PI_RHS(qmc_data)
    global pi_matvec_data
    matvec_data = [b, qmc_data]
    return pi_matvec_data

def SI_MatVec_data(qmc_data):
    """
    MXV_data(qmc_data)
    ------------------
    This function adds the right side of the linear system to Sam's 
    qmc_data structure so I can pass it to the matrix-vector product.
    """
    b   = SI_RHS(qmc_data)
    global si_matvec_data
    matvec_data = [b, qmc_data]
    return si_matvec_data

def SI_MatVec(phi_in):
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
    Nx          = qmc_data.Nx
    G           = qmc_data.G
    Nv          = Nx*G
    phi_s       = np.reshape(phi_in,(Nv,1))
    phi_f       = qmc_data.phi_f
    
    qmc_data.source = np.zeros((Nx,G))
    axv = phi_in - SI_Map(phi_f, phi_s, qmc_data)
    

    return axv


def PI_MatVec(phi_in):
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
    Nx          = qmc_data.Nx
    G           = qmc_data.G
    Nv          = Nx*G
    phi_in      = np.reshape(phi_in,(Nv,1))

    qmc_data.source = np.zeros((Nx,G))
    axv = phi_in - PI_Map(phi_in, qmc_data)
    

    return axv








