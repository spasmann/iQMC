#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon Apr  4 18:07:27 2022

@author: sampasmann

2-Group, Uranium-235, with H2O reflector in Slab geometry. Data taken from
"Analytical Benchmark Test Set For Criticality Code Verification"
"""

import numpy as np

def PUa_H2O_1_0_SL_data(Nx=10):    
    G       = 1
    R       = np.array([1.317862, 2.849725])
    RB      = 2.849725
    LB      = -2.849725
    sigt    = np.empty((Nx,G))
    sigs    = np.empty((Nx,G,G))
    sigf    = np.empty((Nx,G))
    siga    = np.empty((Nx,G))
    chi     = np.empty((Nx,G))
    nu      = np.empty((Nx,G))
    dx      = (RB-LB)/Nx
    xspan   = np.linspace(LB+dx/2,RB-dx/2,num=Nx)
    count   = 0
    for x in xspan:
        if (LB <= x <= -1.317862):
            sigt[count,:]   = 0.32640
            sigs[count,:,:] = 0.293760
            sigf[count,:]   = 0.0
            siga[count,:]   = 0.032640
            chi[count,:]    = 0.0
            nu[count,:]     = 0.0
        elif (-1.317862 < x < 1.317862):
            sigt[count,:]   = 0.32640
            sigs[count,:,:] = 0.225216
            sigf[count,:]   = 0.081600
            siga[count,:]   = 0.081600+0.019584
            chi[count,:]    = 1.0
            nu[count,:]     = 3.24
        elif (1.317862 <= x <= RB):
            sigt[count,:]   = 0.32640
            sigs[count,:,:] = 0.293760
            sigf[count,:]   = 0.0
            siga[count,:]   = 0.032640
            chi[count,:]    = 0.0
            nu[count,:]     = 0.0
        count += 1
        
        
    return sigt, sigs, sigf, siga, chi, nu, G