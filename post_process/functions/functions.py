#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Jun 23 10:46:45 2022

@author: sampasmann
"""

import numpy as np
labeimport h5py
import os
from src.init_files.reeds_solution import reeds_mcdc_sol, reeds_sol
from src.init_files.mg_init import MultiGroupInit

def ReduceFlux(phi, NxRef):
    """
    This function averages over the spatial cells until phi is of the desired 
    length (NxRef)
    
    Parameters
    ----------
    phi : Scalar flux input, the length must be a power of 2 of NxRef
    NxRef : The desired length of the scalar flux

    Returns
    -------
    phi_new : 

    """
    phi_new = np.copy(phi)
    Nx = phi.shape[0]
    I = int(np.log(Nx/NxRef)/np.log(2))
    for i in range(I):
        left_cells = phi_new[0:Nx-1:2,:]
        right_cells = phi_new[1:Nx:2,:]
        phi_new = (right_cells + left_cells)*0.5
        Nx = phi_new.shape[0]
    
    assert (phi_new.shape[0] == NxRef)

    return phi_new

def RelError(phi, sol, order=np.inf):
    RelError = abs(phi - sol)/sol
    return np.linalg.norm(RelError, ord=order)


def AbsError(phi, sol, order=np.inf):
    AbsError = abs(phi-sol)
    return np.linalg.norm(AbsError, ord=order)

def PlotLine(Nvals=np.array((2**10)), Nx=80, generator="halton", problem = "reeds_data", nproc=64):
    #sol = reeds_mcdc_sol()
    data = MultiGroupInit(numGroups=12,Nx=10)
    sol = data.true_flux
    line = np.zeros(len(Nvals))
    path = os.getcwd()+"/../saved_data/"
    NxRef = sol.shape[0]
    G = sol.shape[1]
    
    for i in range(len(Nvals)):
        file = (problem+"-"+generator+"-"+str(Nvals[i])+"-"+str(Nx)+"-"+str(nproc))
        f = h5py.File(path+file, 'r')    
        phi = f['phi_avg'][:]
        f.close()
        phi = np.reshape(phi, (Nx,G))
        phi = ReduceFlux(phi, NxRef)
        error = RelError(phi, sol)
        line[i] = error

    return line
