#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon Apr  4 18:07:27 2022

@author: sampasmann

2-Group, Uranium-235 bare research reactor, Data taken from
"Analytical Benchmark Test Set For Criticality Code Verification"
"""

import numpy as np

def U_2_0_data(Nx=10):    
    G       = 2
    ###########################################
    Nu      = np.array([(2.5,2.7)])
    ###########################################
    Sig_f   = np.array([(0.06912, 0.06192)])
    ###########################################
    Sig_c   = np.array([(0.01344, 0.00384)])
    ###########################################
    Sig_s   = np.array([(0.26304, 0.0720), 
                        (0.00000, 0.078240)])
    ###########################################
    Sig_t   = np.array([(0.3456, 0.2160)])
    ###########################################
    X       = np.array([(0.425, 0.575)])
    ###########################################
    Sig_a   = Sig_c + Sig_f
    
    sigt    = np.tile(Sig_t, (Nx,1))
    sigs    = np.tile(Sig_s, (Nx,1,1))
    sigf    = np.tile(Sig_f, (Nx,1))
    siga    = np.tile(Sig_a, (Nx,1))
    chi     = np.tile(X, (Nx,1))
    nu      = np.tile(Nu, (Nx,1))
    
    return sigt, sigs, sigf, siga, chi, nu, G