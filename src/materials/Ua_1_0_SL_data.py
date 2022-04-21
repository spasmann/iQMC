#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon Apr  4 17:59:34 2022

@author: sampasmann
"""

import numpy as np

def Ua_1_0_SL_data(Nx=100):

    G = 1
    dtype = np.float64
    sigt = np.array((0.32640),dtype)
    sigs = np.array((0.248064),dtype)
    sigc = np.array((0.013056),dtype)
    sigf = np.array((0.065280),dtype)
    siga = sigf + sigc
    nu = np.array((2.70),dtype)
    chi = np.array((1.0),dtype)
    R = np.array([2.872934])
    
    sigt = np.tile(sigt,(Nx,G))
    sigs = np.tile(sigs,(Nx,G,G))
    sigc = np.tile(sigc,(Nx,G))
    sigf = np.tile(sigf,(Nx,G))
    siga = np.tile(siga, (Nx,G))
    nu = np.tile(nu,(Nx,G))
    chi = np.tile(chi,(Nx,G))
    
    return sigt, sigs, siga, sigf, chi, nu, G