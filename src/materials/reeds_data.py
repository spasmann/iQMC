#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon Mar 28 13:22:41 2022

@author: sampasmann
"""
import numpy as np
from time import process_time
import numba

@numba.njit
def reeds_data(Nx=1000):
    G = 1 # number of energy groups
    LB = -8
    RB = 8
    sigt = np.empty((Nx,G))
    sigs = np.empty((Nx,G))
    source = np.empty((Nx,G))
    xspan = range(LB,RB,Nx)
    count = 0
    for x in xspan:
        if (x < -6):
            sigt[count] = 1.0
            sigs[count] = 0.9
            source[count] = 0.0
        elif (-6 <= x < -5):
            sigt[count] = 1.0
            sigs[count] = 0.9
            source[count] = 1.0
        elif (-5 <= x < -3):
            sigt[count] = 0.0
            sigs[count] = 0.0
            source[count] = 0.0
        elif (-3 <= x < -2):
            sigt[count] = 5.0
            sigs[count] = 0.0
            source[count] = 0.0
        elif (-2 <= x < 2):
            sigt[count] = 50.0
            sigs[count] = 0.0
            source[count] = 50.0
        elif (2 <= x < 3):
            sigt[count] = 5.0
            sigs[count] = 0.0
            source[count] = 0.0
        elif (3<= x < 5):
            sigt[count] = 0.0
            sigs[count] = 0.0
            source[count] = 0.0
        elif (5<= x < 6):
            sigt[count] = 1.0
            sigs[count] = 0.9
            source[count] = 1.0
        elif (6<= x):
            sigt[count] = 1.0
            sigs[count] = 0.9
            source[count] = 0.0
        count += 1
        
    siga = (sigt - sigs)
    
    return sigt, sigs, siga, source, G

if __name__ == "__main__":
    start = process_time()
    sigt, sigs, siga, source, G = reeds_data(1)
    stop = process_time()
    time1 = stop-start
    print("Elapsed time with Compliation:", time1)
    
    start = process_time()
    sigt, sigs, siga, source, G = reeds_data(10000)
    stop = process_time()
    time2 = stop-start
    print("Elapsed time After Compliation:", time2)
    
    print("A {}x speed up".format(round(time1/time2)))