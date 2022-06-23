#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Jun 23 10:46:45 2022

@author: sampasmann
"""

import numpy as np

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
    Nx = phi.size[0]
    I = np.log(Nx/NxRef)/np.log(2)
    for i in range(I-1):
        left_cells = phi[0:2:Nx-1,:]
        right_cells = phi[1:2:Nx,:]
        phi_new = (right_cells + left_cells)*0.5
        Nx = phi.size[0]
        
    return phi_new

def RelError(phi, sol, ord=np.inf):
    RelError = abs(phi - sol)/sol
    return np.linalg.norm(RelError, ord=ord)

def AbsError(phi, sol, ord=np.inf):
    AbsError = abs(phi-sol)
    return np.linalg.norm(AbsError, ord=ord)