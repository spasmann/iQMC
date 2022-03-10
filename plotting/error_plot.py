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
fname1 = "12-random-4096-20"
fname2 = "12-sobol-4096-20"
fname3 = "12-halton-4096-20"

f1 = h5py.File(path+fname1, 'r')
f2 = h5py.File(path+fname2, 'r')
f3 = h5py.File(path+fname3, 'r')

print("Keys: ", list(f1.keys()))

err1 = f1['error'][:]
err2 = f2['error'][:]
err3 = f3['error'][:]
itt1 = f1['itt'][...]
itt2 = f2['itt'][...]
itt3 = f3['itt'][...]

ylabel = r'$||\frac{\phi_i - \phi}{\phi}||_\infty$'

plt.figure(dpi=200, figsize=(8,5))
plt.plot(range(itt1), err1[:], label='random')
plt.plot(range(itt3), err2[:], label='halton')
plt.plot(range(itt2), err3[:], label='sobol')
plt.grid()
plt.yscale('log')
plt.xlabel('Iteration')
plt.ylabel(ylabel)
plt.xticks(range(0,itt1+1,2))
plt.legend()

"""
# difference in norms between iterations
plt.figure(dpi=200)
plt.plot(range(itt1-1), diff1[:], label='random')
plt.plot(range(itt2-1), diff2[:], label='sobol')
plt.legend()
"""