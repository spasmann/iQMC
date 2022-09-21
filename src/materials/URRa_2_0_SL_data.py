#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon Apr  4 18:07:27 2022

@author: sampasmann

2-Group, Uranium-235 bare research reactor, Data taken from
"Analytical Benchmark Test Set For Criticality Code Verification"
"""

import numpy as np

def URRa_2_0_SL_data(Nx=10):    
    G       = 2
    Nu      = np.array([(2.5,2.5)])
    Sig_f   = np.array([(0.0010484, 0.050632)])
    Sig_c   = np.array([(0.0010046, 0.025788)])
    
    Sig_s   = np.array([(0.0, 0.029227),
                        (2.44383, 0.62568)])
    Sig_s   = np.flip(Sig_s,axis=1)
    
    Sig_t   = np.array([(0.65696, 2.52025)])
    Sig_a   = Sig_c + Sig_f
    X       = np.array([(1.0,0.0)])
    R       = np.array([7.566853])
    
    sigt    = np.tile(Sig_t, (Nx,1))
    sigs    = np.tile(Sig_s, (Nx,1,1))
    sigf    = np.tile(Sig_f, (Nx,1))
    siga    = np.tile(Sig_a, (Nx,1))
    chi     = np.tile(X, (Nx,1))
    nu      = np.tile(Nu, (Nx,1))
    
        
    return sigt, sigs, sigf, siga, chi, nu, G