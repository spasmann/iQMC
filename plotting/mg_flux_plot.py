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
fname2 = "12-sobol-1024-20"

f2 = h5py.File(path+fname2, 'r')

print("Keys: ", list(f2.keys()))

phi2 = f2['phi_avg'][:]
true2 = f2['true_flux'][:]
RB = f2['RB'][...]
Nx = f2['Nx'][...]
midpoints = np.linspace(0,RB,Nx)

plt.figure(dpi=200,figsize=(8,5))
for G in range(phi2.shape[1]):
    plt.plot(midpoints, phi2[:][:,G], label=G)
    plt.plot(midpoints, true2[:][:,G],'k--')
plt.xlabel('Midoints')
plt.ylabel('Cell Averaged Scalar Flux')
#plt.legend()
