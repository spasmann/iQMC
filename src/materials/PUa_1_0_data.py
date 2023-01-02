#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon Apr  4 17:59:34 2022

@author: sampasmann
"""

import numpy as np

def PUa_1_0_data(Nx=20):

    G = 1
    dtype = np.float64
    
    sigt = np.array((0.32640),dtype)
    sigs = np.array((0.225216),dtype)
    sigc = np.array((0.019584),dtype)
    sigf = np.array((0.081600),dtype)
    siga = sigf + sigc
    nu = np.array((2.84),dtype)
    chi = np.array((1.0),dtype)
    
    sigt = np.tile(sigt,(Nx,G))
    sigs = np.tile(sigs,(Nx,G,G))
    sigc = np.tile(sigc,(Nx,G))
    sigf = np.tile(sigf,(Nx,G))
    siga = np.tile(siga, (Nx,G))
    nu = np.tile(nu,(Nx,G))
    chi = np.tile(chi,(Nx,G,G))
    
    return sigt, sigs, siga, sigf, chi, nu, G