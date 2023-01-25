#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Jun 23 10:46:45 2022

@author: sampasmann
"""

import sys, os
sys.path.append(os.getcwd()+"/../../")
import numpy as np
import h5py
import os
from src.input_files.reeds_solution import reeds_mcdc_sol, reeds_sol
from src.input_files.mg_init import MultiGroupInit

# =============================================================================
# Spatial Averaging Technique
# =============================================================================
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


# =============================================================================
# Garcia Angular Flux Retrieval
# =============================================================================
def garcia_angle_bins():
    angles = np.array((-0.05))
    angles = np.append(angles, np.arange(-0.1,-1.1,step=-0.1))
    angles = np.append(angles, 0.05)
    angles = np.append(angles, np.arange(0.1,1.1,step=0.1))
    return angles

def garcia_angular_flux_sol():
    sol_left=np.array((8.97797e-01, 8.87836e-01, 8.69581e-01, 8.52299e-01, 
                       8.35503e-01, 8.18996e-01, 8.02676e-01, 7.86493e-01, 
                       7.70429e-01, 7.54496e-01, 7.38721e-01))
    sol_right= np.array((1.02202e-01, 1.12164e-01, 1.30419e-01, 1.47701e-01, 
                         1.64497e-01, 1.81004e-01, 1.97324e-01, 2.13507e-01, 
                         2.29571e-01, 2.45504e-01, 2.61279e-01))
    return (sol_left, sol_right)

def SN_Sweep(angles, qmc_data):
    #
    # data from angles
    #
    Na2 = angles.size
    Na  = int(Na2*0.5)
    #
    # data from iQMC
    #
    Nx              = qmc_data.Nx
    dx              = qmc_data.mesh.dx
    phi             = qmc_data.tallies.phi_avg
    fixed_source    = qmc_data.fixed_source
    sigt            = qmc_data.material.sigt
    sigs            = qmc_data.material.sigs
    psi             = np.zeros((Na2,Nx))
    psi[:,0]        = qmc_data.phi_left
    psi[:,-1]       = qmc_data.phi_right
    #
    # data for SN sweep
    #
    source_total    = (0.5 * sigs[:,0,0] * phi[:,0] + fixed_source[:,0])
    forward_angles  = angles[Na:Na2]
    # backward_angles = angles[:Na]
    #
    # These next terms should technically be in the for loop to index the XS
    # also we only use the forward angles because the formula calls for 
    # the abs(\mu) and the forward/backward angles are the same in this case
    # SN-DD Denominator
    vfl             = (forward_angles / dx) + sigt[0,0]*0.5
    vfl             = 1.0 / vfl
    # SN-DD part of numerator
    vfr             = (forward_angles / dx) - sigt[0,0]*0.5
    #
    # forward sweep
    #
    for i in range(1,Nx):
        psi[Na:Na2, i]  = np.copy(psi[Na:Na2, i-1])
        psi[Na:Na2, i] *= vfr
        psi[Na:Na2, i] += source_total[i]
        psi[Na:Na2, i] *= vfl
    #
    # backward sweep
    #
    for i in range(Nx-2,-1,-1):
        psi[:Na, i]  = np.copy(psi[:Na, i+1])
        psi[:Na, i] *= vfr
        psi[:Na, i] += source_total[i]
        psi[:Na, i] *= vfl
    
    return psi


def MOC_Sweep(angles, qmc_data):
    #
    # data from angles
    #
    Na2             = angles.size
    Na              = int(Na2*0.5)
    forward_angles  = angles[Na:Na2]
    backward_angles = angles[:Na]
    #
    # data from iQMC
    #
    Nx              = qmc_data.Nx
    dx              = qmc_data.mesh.dx
    phi             = qmc_data.tallies.phi_avg
    fixed_source    = qmc_data.fixed_source
    xspan           = qmc_data.mesh.midpoints
    sigt            = qmc_data.material.sigt
    sigs            = qmc_data.material.sigs
    psi             = np.zeros((Na2,Nx))
    psi[:,0]        = qmc_data.phi_left
    psi[:,-1]       = qmc_data.phi_right
    #
    # Forward sweep
    #
    for i in range(1,Nx):
        psi[Na:Na2, i]  = np.copy(psi[Na:Na2, i-1])
        source_total    =  (0.5 * sigs[i,0,0] * phi[i,0] + fixed_source[i,0])
        psi[Na:Na2, i] -= source_total / sigt[i,0]
        psi[Na:Na2, i] *= np.exp(-sigt[i,0] * dx / forward_angles)
        psi[Na:Na2, i] += source_total / sigt[i,0]
    #
    # backward sweep
    #
    for i in range(Nx-2,-1,-1):
        psi[:Na, i]  = np.copy(psi[Na:Na2, i+1])
        source_total =  (0.5 * sigs[i,0,0] * phi[i,0] + fixed_source[i,0])
        psi[:Na, i] -= source_total / sigt[i,0]
        psi[:Na, i] *= np.exp(-sigt[i,0] * dx / forward_angles)
        psi[:Na, i] += source_total / sigt[i,0]
    
    return psi

# =============================================================================
# Misc
# =============================================================================
def RelError(phi, sol, order=np.inf):
    RelError = abs(phi - sol)/sol
    return np.linalg.norm(RelError, ord=order)


def AbsError(phi, sol, order=np.inf):
    AbsError = abs(phi-sol)
    return np.linalg.norm(AbsError, ord=order)

def PlotLine(Nvals=np.array((2**10)), Nx=80, generator="halton", problem = 
             "reeds_data", nproc=64):
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

if (__name__ == "__main__"):
    from src.input_files.garcia_init import GarciaInit
    data = GarciaInit(N=2**6, Nx=4)
    angles = garcia_angle_bins()
    SN_Sweep(angles, data)