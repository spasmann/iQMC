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
fname1 = "reeds_data-halton-16384-180"
fname2 = "reeds_data-halton-4096-180"
fname3 = "reeds_data-halton-1024-180"
<<<<<<< HEAD
fname5 = "12-halton-1024-20"
=======
fname5 = "reeds_data-halton-2000-400"
>>>>>>> 493a1141ed56ef9e67784f5c8e613194ca5045cc
files = [fname5]

plt.figure(dpi=200)
count = 0
for file in files:
    f = h5py.File(path+file, 'r')
    print("Keys: ", list(f.keys()))
    phi = f['phi_avg'][:]
    true = f['true_flux'][:]
    RB = f['RB'][...]
    LB = f['LB'][...]
    Nx = f['Nx'][...]
    N = f['N'][...]
    generator = f['generator'][...]
    midpoints = np.linspace(LB,RB,Nx)
    for G in range(phi.shape[1]):
        plt.plot(midpoints, phi[:][:,G], label=files[count])

    f.close()
    count += 1
plt.plot(midpoints, true[:],'--', label="true")  

plt.xlabel(r'Spatial Position $x$')
plt.ylabel(r'Cell Averaged Scalar Flux $\phi$')
#plt.legend()
