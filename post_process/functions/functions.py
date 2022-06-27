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
    Nx = len(phi)
    I = int(np.log(Nx/NxRef)/np.log(2))
    for i in range(I):
        left_cells = phi_new[0:Nx-1:2,:]
        right_cells = phi_new[1:Nx:2,:]
        phi_new = (right_cells + left_cells)*0.5
        Nx = len(phi)
        
    return phi_new

def RelError(phi, sol, order=np.inf):
    assert (phi.shape == sol.shape)
    RelError = abs(phi - sol)/sol
    #return np.linalg.norm(RelError)
    return RelError.max()


def AbsError(phi, sol, order=np.inf):
    assert (phi.shape == sol.shape)
    AbsError = abs(phi-sol)
    #return np.linalg.norm(AbsError)
    return AbsError.max()
