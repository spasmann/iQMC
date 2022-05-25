#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue May  3 14:50:16 2022

@author: sampasmann
"""

import matplotlib.pyplot as plt
import matplotlib
import numpy as np
import h5py

path = "../saved_data/"
fname1 = "reeds_data-latin_hypercube-2048-320"
fname2 = "reeds_data-halton-32768-1280"
fname3 = "reeds_data-halton-2000-400"
fname4 = "reeds_data-halton-2000-560"

f1 = h5py.File(path+fname1, 'r')
f2 = h5py.File(path+fname4, 'r')
nx1 = f1['Nx'][...]
nx2 = f2['Nx'][...]
N1 = f1['N'][...]
N2 = f2['N'][...]
phi1 = f1['phi_avg'][...]
true1 = f1['true_flux'][...]
phi2 = f2['phi_avg'][...]
true2 = f2['true_flux'][...]

error1 = f1['error'][...]
RE1 = abs((phi1-true1))
RE2 = abs((phi2-true2))

LB = -8
RB = 8
dx1 = (RB-LB)/nx1
xspan1 = np.linspace(LB+dx1/2, RB-dx1/2, nx1)
dx2 = (RB-LB)/nx2
xspan2 = np.linspace(LB+dx2/2, RB-dx2/2, nx2)

plt.figure()
plt.plot(xspan1,RE1,label='N={},Nx={}'.format(N1,nx1))
#plt.plot(xspan2,RE2,label='N={},Nx={}'.format(N2,nx2))
plt.title("Relative Error")
plt.xlabel("Spatial Domain".format(nx1))
plt.legend()

f1.close()
f2.close()