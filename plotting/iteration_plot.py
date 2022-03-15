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

phi1 = f1['phi_avg'][:]
phi2 = f2['phi_avg'][:]
phi3 = f3['phi_evg'][:]
itt1 = f1['itt'][...]
itt2 = f2['itt'][...]
itt3 = f3['itt'][...]

err1 = np.linalg.norm(norm1[:itt1-1] - norm1[1:])
err2 = norm2[:itt2-1] - norm2[1:]
err2 = norm3[:itt3-1] - norm3[1:]

plt.figure(dpi=200, figsize=(8,5))
plt.plot(range(itt1), norm1[:], label='random')
plt.plot(range(itt3), norm3[:], label='halton')
plt.plot(range(itt2), norm2[:], label='sobol')
plt.grid()
plt.yscale('log')
plt.xticks(range(0,itt1+1,2))
plt.legend()

"""
# difference in norms between iterations
plt.figure(dpi=200)
plt.plot(range(itt1-1), diff1[:], label='random')
plt.plot(range(itt2-1), diff2[:], label='sobol')
plt.legend()
"""