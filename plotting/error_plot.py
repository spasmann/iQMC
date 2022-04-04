#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon Mar  7 14:36:21 2022

@author: sampasmann
"""

import matplotlib.pyplot as plt
import matplotlib
import numpy as np
import h5py
matplotlib.rcParams.update({'font.size': 16})

path = "../saved_data/"
fname1 = "reeds_data-halton-1024-180"

files = [fname1]

plt.figure(dpi=200)
for file in files:
    f = h5py.File(path+file, 'r')
    print("Keys: ", list(f.keys()))
    err = f['error'][:]
    itt = f['itt'][...]
    phi = f['phi_avg'][:]
    N = f['N'][...]
    generator = f['generator'][...]
    for G in range(phi.shape[1]):
        plt.plot(range(itt), err[:], label=generator)
    f.close()    

ylabel = r'$||\frac{\phi_i - \phi}{\phi}||_\infty$'
plt.grid()
plt.yscale('log')
plt.xlabel('Iteration')
plt.ylabel(ylabel)
plt.xticks(range(0,itt+1,2))
plt.legend()


