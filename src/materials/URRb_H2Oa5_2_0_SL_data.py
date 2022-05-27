#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon Apr  4 18:07:27 2022

@author: sampasmann

2-Group, Uranium-235, with H2O reflector in Slab geometry. Data taken from
"Analytical Benchmark Test Set For Criticality Code Verification"
"""

import numpy as np

def URRb_H2Oa5_2_0_SL_data(Nx=10):    
    G = 2

    Nu = np.array([(2.5,2.5),(0.0,0.0)])
    Sig_f = np.array([(0.000836, 0.029564),(0.0,0.0)])
    Sig_c = np.array([(0.001104, 0.024069),(0.00074, 0.018564)])
    Sig_s = np.array([((2.9183,0.000767),(0.04635,0.83892)),
                       ((2.9676,0.000336),(0.04749,0.83975))])
    Sig_t = np.array([(0.88721, 2.9727),(0.88798, 2.9865)])
    Sig_a = Sig_c + Sig_f
    X = np.array([(1.0,0.0),(0.0,0.0)])
    R = np.array([6.696802, 7.822954])
    
    #Sig_s = np.flip(Sig_s,1)
    
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
            sigt[count,:] = Sig_t[1,:]
            sigs[count,:,:] = Sig_s[1,:]
            sigf[count,:] = Sig_f[1,:]
            siga[count,:] = Sig_a[1,:]
            chi[count,:]  = X[1,:]
            nu[count,:]   = Nu[1,:]
        elif (-R[0] < x < R[0]):
            sigt[count,:] = Sig_t[0,:]
            sigs[count,:,:] = Sig_s[0,:]
            sigf[count,:] = Sig_f[0,:]
            siga[count,:] = Sig_a[0,:]
            chi[count,:]  = X[0,:]
            nu[count,:]   = Nu[0,:]
        elif ( R[0] <= x <= R[1]):
            sigt[count,:] = Sig_t[1,:]
            sigs[count,:,:] = Sig_s[1,:]
            sigf[count,:] = Sig_f[1,:]
            siga[count,:] = Sig_a[1,:]
            chi[count,:]  = X[1,:]
            nu[count,:]   = Nu[1,:]
        count += 1
        
        
    return sigt, sigs, sigf, siga, chi, nu, G