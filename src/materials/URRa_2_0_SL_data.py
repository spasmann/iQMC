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
    G = 2
    Nu = np.array([(2.5,2.5)])
    Sig_f = np.array([(0.050632,0.0010484)])
    Sig_c = np.array([(0.025788,0.0010046)])
    Sig_s = np.array([(2.44383,0.029227),(0.0, 0.62568)])
    Sig_t = np.array([(2.52025,0.65696)])
    Sig_a = Sig_c + Sig_f
    X = np.array([(0.0,1.0)])
    R = np.array([7.566853])
    
    sigt = np.zeros((Nx,G))
    #sigs = np.zeros((Nx,G,G))
    sigf = np.zeros((Nx,G))
    siga = np.zeros((Nx,G))
    chi = np.zeros((Nx,G))
    nu = np.zeros((Nx,G))
    xspan = np.linspace(-R[0],R[0],num=Nx)
    count = 0
    for x in xspan:
        if (-R[0] <= x <= R[0]):
            sigt[count,:] = Sig_t[:]
            #sigs[count,:,:] = Sig_s[:]
            sigf[count,:] = Sig_f[:]
            #siga[count,:] = Sig_a[:]
            chi[count,:]  = X[:]
            nu[count,:]   = Nu[:]
        count += 1

        
    return sigt, Sig_s, sigf, siga, chi, nu, G