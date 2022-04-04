#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon Mar 28 13:22:41 2022

@author: sampasmann
"""
import numpy as np

def reeds_data(Nx=1000):
    G = 1 # number of energy groups
    sigt = np.zeros((Nx,1))
    sigs = np.zeros((Nx,1))
    source = np.zeros((Nx,1))
    xspan = np.linspace(-8.0,8.0,num=Nx)
    count = 0
    for x in xspan:
        if (x < -6.0):
            sigt[count] = 1.0
            sigs[count] = 0.9
            source[count] = 0.0
        elif (-6.0 <= x < -5.0):
            sigt[count] = 1.0
            sigs[count] = 0.9
            source[count] = 1.0
        elif (-5.0 <= x < -3.0):
            sigt[count] = 0.0
            sigs[count] = 0.0
            source[count] = 0.0
        elif (-3.0 <= x < -2.0):
            sigt[count] = 5.0
            sigs[count] = 0.0
            source[count] = 0.0
        elif (-2.0 <= x < 2.0):
            sigt[count] = 50.0
            sigs[count] = 0.0
            source[count] = 50.0
        elif (2.0 <= x < 3.0):
            sigt[count] = 5.0
            sigs[count] = 0.0
            source[count] = 0.0
        elif (3.0<= x < 5.0):
            sigt[count] = 0.0
            sigs[count] = 0.0
            source[count] = 0.0
        elif (5.0<= x < 6.0):
            sigt[count] = 1.0
            sigs[count] = 0.9
            source[count] = 1.0
        elif (6.0<= x):
            sigt[count] = 1.0
            sigs[count] = 0.9
            source[count] = 0.0
        count += 1
        
    siga = (sigt - sigs)
    
    return sigt, sigs, siga, source, G

