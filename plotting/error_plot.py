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


path = "../saved_data/"
fname1 = "reeds_data-halton-2048-180"
fname2 = "reeds_data-sobol-2048-180"
fname3 = "reeds_data-latin_hypercube-2048-180"
fname4 = "reeds_data-random-2048-180"

fname5 = "reeds_data-halton-32768-32"
files = [fname5]


plt.figure(dpi=200,figsize=(15,8))
count = 0
for file in files:
    f = h5py.File(path+file, 'r')
    print("Keys: ", list(f.keys()))
    err = f['error'][:]
    itt = f['itt'][...]
    phi = f['phi_avg'][:]
    N = f['N'][...]
    generator = f['generator'][...]
    plt.plot(range(itt), err[:], label=files[count])
    f.close()    
    count += 1

plt.legend(loc="upper right")
matplotlib.rcParams.update({'font.size': 16})
ylabel = r'$||\frac{\phi_i - \phi}{\phi}||_\infty$'
plt.grid()
plt.yscale('log')
plt.xlabel('Iteration')
plt.ylabel(ylabel)
plt.xticks(range(0,itt+1,2))




