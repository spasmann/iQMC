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
fname1 = "reeds_data-sobol-1024-180"
fname2 = "reeds_data-sobol-2048-180"
fname3 = "reeds_data-halton-2048-180"
files = [fname1, fname2, fname3]

plt.figure(dpi=200)
for file in files:
    f = h5py.File(path+file, 'r')
    print("Keys: ", list(f.keys()))
    phi = f['phi_avg'][:]
    RB = f['RB'][...]
    LB = f['LB'][...]
    Nx = f['Nx'][...]
    N = f['N'][...]
    generator = f['generator'][...]
    midpoints = np.linspace(LB,RB,Nx)
    for G in range(phi.shape[1]):
        plt.plot(midpoints, phi[:][:,G], label=generator)
    f.close()    
    

plt.xlabel(r'Spatial Position $x$')
plt.ylabel(r'Cell Averaged Scalar Flux $\phi$')
plt.legend()
