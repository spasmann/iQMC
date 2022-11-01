#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon Apr  4 17:59:34 2022

@author: sampasmann
"""

import numpy as np
import os

def PUa_1_0_SL_data(Nx=100):

    G = 1
    dtype = np.float64
    
    sigt = np.array((0.32640),dtype)
    sigs = np.array((0.225216),dtype)
    sigc = np.array((0.019584),dtype)
    sigf = np.array((0.081600),dtype)
    siga = sigf + sigc
    nu = np.array((3.24),dtype)
    chi = np.array((1.0),dtype)
    
    sigt = np.tile(sigt,(Nx,G))
    sigs = np.tile(sigs,(Nx,G,G))
    sigc = np.tile(sigc,(Nx,G))
    sigf = np.tile(sigf,(Nx,G))
    siga = np.tile(siga, (Nx,G))
    nu = np.tile(nu,(Nx,G))
    chi = np.tile(chi,(Nx,G,G))
    
    """
    G = 2
    Nu = np.array([(3.24, 3.24)])
    Sig_f = np.array([(0.081600,0.081600)])
    Sig_c = np.array([(0.019584,0.019584)])
    Sig_s = np.array([(0.225216/2,0.225216/2),(0.225216/2, 0.225216/2)])
    Sig_t = np.array([(0.32640,0.32640)])
    Sig_a = Sig_c + Sig_f
    X = np.array([(0.5,0.5)])
    R = np.array([1.853722])
    
    sigt = np.zeros((Nx,G))
    sigs = np.zeros((Nx,G,G))
    sigf = np.zeros((Nx,G))
    siga = np.zeros((Nx,G))
    chi = np.zeros((Nx,G))
    nu = np.zeros((Nx,G))
    for count in range(Nx):
        sigt[count,:] = Sig_t[:]
        sigs[count,:,:] = Sig_s[:]
        sigf[count,:] = Sig_f[:]
        #siga[count,:] = Sig_a[:]
        chi[count,:]  = X[:]
        nu[count,:]   = Nu[:]
    """
    return sigt, sigs, siga, sigf, chi, nu, G