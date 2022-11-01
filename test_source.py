#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue Nov  1 10:50:32 2022

@author: sampasmann
"""

#from src.input_files.mg_init import MultiGroupInit
from src.input_files.URRa_2_0_SL_init import URRa_2_0_SL_init
from src.functions.source import GetSource, scattering_source, fission_source
import numpy as np

def test_source():
    N = 1
    Nx = 1
    generator = "sobol"
    #mg_data = MultiGroupInit(N=N, Nx=Nx, generator=generator)
    data = URRa_2_0_SL_init(N=N, Nx=Nx, generator=generator)
    
    #mg_sigt = mg_data.material.sigt
    #mg_sigs = mg_data.material.sigs
    
    #sigt    = data.material.sigt[0,:]
    sigs    = data.material.sigs[0,:,:]
    sigf    = data.material.sigf[0,:]
    chi     = data.material.chi[0,:]
    nu      = data.material.nu[0,:]
    keff    = data.keff
    
    np.random.seed()
    phi0    = np.random.random((1,2))
    phi0_1  = phi0[0,0]
    phi0_2  = phi0[0,1] 
    
    Sig11 = (sigs[0,0] + chi[0]*nu[0]*sigf[0]/keff) # transfer from 1 to 1
    Sig12 = (sigs[0,1] + chi[0]*nu[1]*sigf[1]/keff) # transfer from 2 to 1 (down scatter)
    Sig21 = (sigs[1,0] + chi[1]*nu[0]*sigf[0]/keff) # transfer from 1 to 2 (up scatter)
    Sig22 = (sigs[1,1] + chi[1]*nu[1]*sigf[1]/keff) # transfer from 2 to 2
    
    scatter = np.array(((sigs[0,0]*phi0_1 + sigs[0,1]*phi0_2),
                        (sigs[1,0]*phi0_1 + sigs[1,1]*phi0_2)))
    fission = np.array(((chi[0]*nu[0]*sigf[0]/keff*phi0_1 + chi[0]*nu[1]*sigf[1]/keff*phi0_2),
                        (chi[1]*nu[0]*sigf[0]/keff*phi0_1 + chi[1]*nu[1]*sigf[1]/keff*phi0_2)))
    
    
    phi1  = np.array(((Sig11*phi0_1 + Sig12*phi0_2),
                      (Sig21*phi0_1 + Sig22*phi0_2)))
    scatter_src = scattering_source(0, phi0, data.material)
    fission_src = fission_source(0, phi0, keff, data.material)
    phi1_src = GetSource(phi0, data,  phi_avg_f=phi0)[0,:]
    
    print(abs(scatter-scatter_src))
    print(abs(fission-fission_src))
    print(abs(phi1-phi1_src))
    
    assert(abs(scatter-scatter_src).all() == 0.0)
    assert(abs(fission-fission_src).all() == 0.0)
    assert(abs(phi1-phi1_src).all() == 0.0)

if (__name__ == "__main__"):
    test_source()



