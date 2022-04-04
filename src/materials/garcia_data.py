#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon Apr  4 14:44:55 2022

@author: sampasmann
"""

import numpy as np

def garcia_data(mesh, Nx=1000):
    G = 1
    sigt = np.ones((Nx, G))
    c = 1.0
    sigs = np.exp(-mesh.midpoints/c)
    sigs = np.reshape(sigs, (Nx, G))
    siga = sigt - sigs
    
    return sigt, sigs, siga, G