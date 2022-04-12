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
fname1 = "reeds_data-random-4096-180"
fname2 = "moment_match_random"
fname3 = "reeds_data-halton-4096-180"
fname4 = "moment_match_halton"
files = [fname1, fname2, fname3, fname4]


plt.figure(dpi=200)
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

plt.legend()
matplotlib.rcParams.update({'font.size': 16})
ylabel = r'$||\frac{\phi_i - \phi}{\phi}||_\infty$'
plt.grid()
plt.yscale('log')
plt.xlabel('Iteration')
plt.ylabel(ylabel)
plt.xticks(range(0,itt+1,2))




