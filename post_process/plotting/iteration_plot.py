#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon Mar  7 14:36:21 2022

@author: sampasmann
"""

import matplotlib.pyplot as plt
import numpy as np
import h5py

path = "../saved_data/"
fname1 = "12-random-1024-20"
fname2 = "12-sobol-1024-20"
fname3 = "12-halton-1024-20"

f1 = h5py.File(path+fname1, 'r')
f2 = h5py.File(path+fname2, 'r')
f3 = h5py.File(path+fname3, 'r')

print("Keys: ", list(f1.keys()))

dphi1 = f1['delta_flux'][:]
dphi2 = f2['delta_flux'][:]
dphi3 = f3['delta_flux'][:]
itt1 = f1['itt'][...]
itt2 = f2['itt'][...]
itt3 = f3['itt'][...]

plt.figure(dpi=200, figsize=(8,5))
plt.plot(range(itt1), dphi1[:], label='random')
plt.plot(range(itt2), dphi2[:], label='sobol')
plt.plot(range(itt3), dphi3[:], label='halton')
plt.grid()
plt.yscale('log')
plt.xticks(range(0,itt3+1,2))
plt.legend()
plt.xlabel('Number of Iterations')
plt.ylabel('Relative Residual')

"""
# difference in norms between iterations
plt.figure(dpi=200)
plt.plot(range(itt1-1), diff1[:], label='random')
plt.plot(range(itt2-1), diff2[:], label='sobol')
plt.legend()
"""