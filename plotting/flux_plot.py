#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon Mar  7 14:36:08 2022

@author: sampasmann
"""

import matplotlib.pyplot as plt
import numpy as np
import h5py


path = "../saved_data/"
fname1 = "test_data-random-1024-20"
fname2 = "test_data-sobol-1024-20"

f1 = h5py.File(path+fname1, 'r')
f2 = h5py.File(path+fname2, 'r')
print("Keys: ", list(f1.keys()))

phi1 = f1['phi_avg'][:]
phi2 = f2['phi_avg'][:]
RB = f1['RB'][...]
Nx = f1['Nx'][...]
midpoints = np.linspace(0,RB,Nx)

plt.figure(dpi=200)
for G in range(phi1.shape[1]):
    plt.plot(midpoints, phi1[:][:,G], label='random')
    plt.plot(midpoints, phi2[:][:,G], label='sobol')
plt.xlabel('Midoints')
plt.ylabel('Cell Averaged Scalar Flux')
plt.legend()
