#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon Apr  4 17:59:34 2022

@author: sampasmann
"""

import numpy as np
import os

def pu239_data(Nx=100):
    G = 1
    sigt = np.array((0.32640))
    sigs = np.array((0.225216))
    sigc = np.array((0.019584))
    sigf = np.array((0.081600))
    siga = sigf + sigc
    nu = np.array((3.24))
    
    sigt = np.tile(sigt,(Nx,G))
    sigs = np.tile(sigs,(Nx,G))
    sigc = np.tile(sigc,(Nx,G))
    sigf = np.tile(sigf,(Nx,G))
    siga = np.tile(siga, (Nx,G))
    nu = np.tile(nu,(Nx,G))
    
    
    return sigt, sigs, siga, sigf, nu, G