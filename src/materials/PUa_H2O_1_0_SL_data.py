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
    G = 1
    """
    Nu = np.array([(3.24),(0.0)])
    Sig_f = np.array([(0.081600),(0.0)])
    Sig_c = np.array([(0.019584),(0.032640)])
    Sig_s = np.array([(0.225216), (0.293760)])
    Sig_t = np.array([(0.32640),(0.32640)])
    Sig_a = Sig_c + Sig_f
    X = np.array([(1.0),(0.0)])
    R = np.array([1.317862, 2.849725])

    sigt = np.zeros((Nx,G))
    sigs = np.zeros((Nx,G,G))
    sigf = np.zeros((Nx,G))
    siga = np.zeros((Nx,G))
    chi = np.zeros((Nx,G))
    nu = np.zeros((Nx,G))
    
    xspan = np.linspace(-R[-1],R[-1],num=Nx)
    count = 0
    for x in xspan:
        if (-R[1] <= x <= -R[0]):
            sigt[count,:] = Sig_t[1]
            sigs[count,:,:] = Sig_s[1]
            sigf[count,:] = Sig_f[1]
            siga[count,:] = Sig_a[1]
            chi[count,:]  = X[1]
            nu[count,:]   = Nu[1]
        elif (-R[0] < x < R[0]):
            sigt[count,:] = Sig_t[0]
            sigs[count,:,:] = Sig_s[0]
            sigf[count,:] = Sig_f[0]
            siga[count,:] = Sig_a[0]
            chi[count,:]  = X[0]
            nu[count,:]   = Nu[0]
        elif ( R[0] <= x <= R[1]):
            sigt[count,:] = Sig_t[1]
            sigs[count,:,:] = Sig_s[1]
            sigf[count,:] = Sig_f[1]
            siga[count,:] = Sig_a[1]
            chi[count,:]  = X[1]
            nu[count,:]   = Nu[1]
        count += 1
    """
    R = np.array([1.317862, 2.849725])

    sigt = np.zeros((Nx,G))
    sigs = np.zeros((Nx,G,G))
    sigf = np.zeros((Nx,G))
    siga = np.zeros((Nx,G))
    chi = np.zeros((Nx,G))
    nu = np.zeros((Nx,G))
    
    xspan = np.linspace(-R[-1],R[-1],num=Nx)
    count = 0
    for x in xspan:
        if (-R[1] <= x <= -R[0]):
            sigt[count,:] = 0.32640
            sigs[count,:,:] = 0.293760
            sigf[count,:] = 0.0
            siga[count,:] = 0.032640
            chi[count,:]  = 0.0
            nu[count,:]   = 0.0
            
        elif (-R[0] < x < R[0]):
            sigt[count,:] = 0.32640
            sigs[count,:,:] = 0.225216
            sigf[count,:] = 0.081600
            siga[count,:] = 0.081600+0.019584
            chi[count,:]  = 1.0
            nu[count,:]   = 3.24
            
        elif ( R[0] <= x <= R[1]):
            sigt[count,:] = 0.32640
            sigs[count,:,:] = 0.293760
            sigf[count,:] = 0.0
            siga[count,:] = 0.032640
            chi[count,:]  = 0.0
            nu[count,:]   = 0.0
        count += 1
        
        
    return sigt, sigs, sigf, siga, chi, nu, G