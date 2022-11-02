#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue Nov  1 10:50:32 2022

@author: sampasmann
"""
import sys
sys.path.append("../../")
from src.input_files.U_2_0_SL_init import U_2_0_SL_init
from src.functions.source import GetSource, scattering_source, fission_source
import numpy as np

def test_source():
    # generate test data
    N           = 1
    Nx          = 1
    generator   = "sobol"
    data        = U_2_0_SL_init(N=N, Nx=Nx, generator=generator)
    # pull cross sections
    sigs    = data.material.sigs[0,:,:]
    sigf    = data.material.sigf[0,:]
    chi     = data.material.chi[0,:]
    nu      = data.material.nu[0,:]
    keff    = data.keff
    # set random seed for random flux values
    np.random.seed()
    phi0    = np.random.random((1,2))
    phi0_1  = phi0[0,0]
    phi0_2  = phi0[0,1] 
    # calculate the group transfer rates
    Sig11 = (sigs[0,0] + chi[0]*nu[0]*sigf[0]/keff) # transfer from 1 to 1
    Sig12 = (sigs[0,1] + chi[0]*nu[1]*sigf[1]/keff) # transfer from 2 to 1 (down scatter)
    Sig21 = (sigs[1,0] + chi[1]*nu[0]*sigf[0]/keff) # transfer from 1 to 2 (up scatter)
    Sig22 = (sigs[1,1] + chi[1]*nu[1]*sigf[1]/keff) # transfer from 2 to 2
    phi1  = np.array(((Sig11*phi0_1 + Sig12*phi0_2),# new scalar flux
                      (Sig21*phi0_1 + Sig22*phi0_2))) 
    # calculate the scattering and fission sources separately
    scatter = np.array(((sigs[0,0]*phi0_1 + sigs[0,1]*phi0_2),
                        (sigs[1,0]*phi0_1 + sigs[1,1]*phi0_2)))
    fission = np.array(((chi[0]*nu[0]*sigf[0]/keff*phi0_1 + chi[0]*nu[1]*sigf[1]/keff*phi0_2),
                        (chi[1]*nu[0]*sigf[0]/keff*phi0_1 + chi[1]*nu[1]*sigf[1]/keff*phi0_2)))
    
    # calculate sources using src code functions
    scatter_src = scattering_source(0, phi0, data.material)
    fission_src = fission_source(0, phi0, keff, data.material)
    phi1_src    = GetSource(phi0, data,  phi_avg_f=phi0)[0,:]
    
    # the total, fission, and scattering sources should all match
    #print(abs(scatter-scatter_src))
    #print(abs(fission-fission_src))
    #print(abs(phi1-phi1_src))
    
    assert(np.isclose(abs(scatter-scatter_src).sum(), 0.0))
    assert(np.isclose(abs(fission-fission_src).sum(), 0.0))
    assert(np.isclose(abs(phi1-phi1_src).sum(), 0.0))

if (__name__ == "__main__"):
    test_source()



