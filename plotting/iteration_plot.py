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
fname1 = "garcia_data-random-1024-20"
fname2 = "garcia_data-sobol-1024-20"

f1 = h5py.File(path+fname1, 'r')
f2 = h5py.File(path+fname2, 'r')

print("Keys: ", list(f1.keys()))

norm1 = f1['delta_flux'][:]
norm2 = f2['delta_flux'][:]
itt1 = f1['itt'][...]
itt2 = f2['itt'][...]

diff1 = norm1[:itt1-1] - norm1[1:]
diff2 = norm2[:itt2-1] - norm2[1:]


plt.figure(dpi=200)
plt.plot(range(itt1), norm1[:], label='random')
plt.plot(range(itt2), norm2[:], label='sobol')

plt.figure(dpi=200)
plt.plot(range(itt1-1), diff1[:], label='random')
plt.plot(range(itt2-1), diff2[:], label='sobol')
