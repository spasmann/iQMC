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
fname1 = "reeds_data-halton-32768-80"
fname2 = "reeds_data-halton-32768-160"
fname3 = "reeds_data-halton-32768-320"
fname4 = "reeds_data-halton-32768-640"
fname5 = "reeds_data-halton-32768-1280"
fname6 = "reeds_data-halton-32768-2560"
fname7 = "reeds_data-halton-65536-80"
fname8 = "reeds_data-halton-65536-160"
fname9 = "reeds_data-halton-65536-320"
fname10 = "reeds_data-halton-65536-640"
fname11 = "reeds_data-halton-65536-1280"
fname12 = "reeds_data-halton-65536-2560"
fname13 = "reeds_data-halton-131072-80"
fname14 = "reeds_data-halton-131072-160"
fname15 = "reeds_data-halton-131072-320"
fname16 = "reeds_data-halton-131072-640"
fname17 = "reeds_data-halton-131072-1280"
fname18 = "reeds_data-halton-131072-2560"
fname19 = "reeds_data-halton-262144-80"
fname20 = "reeds_data-halton-262144-160"
fname21 = "reeds_data-halton-262144-320"
fname22 = "reeds_data-halton-262144-640"
fname23 = "reeds_data-halton-262144-1280"

"""
files = [fname1,fname2,fname3,fname4,fname5,fname6,fname7,
         fname8,fname9,fname10,fname11,fname12,fname13,fname14,
         fname15,fname16,fname17,fname18,fname19,fname20,fname21,
         fname22,fname23]
"""
N32768 = [fname1,fname2,fname3,fname4,fname5,fname6]
N65536 = [fname7,fname8,fname9,fname10,fname11,fname12]
N131072 = [fname13,fname14,fname15,fname16,fname17,fname18]
N262144 = [fname19,fname20,fname21,fname22,fname2]


files = [N32768, N65536, N131072, N262144]
numN = len(files)
numNx = 5
N = [32768, 65536, 131072, 262144]
Nx = [80, 80*2, 80*4, 80*8, 80*16]
plt.figure(dpi=200,figsize=(8,5))

count = 0
for i in range(numN):
    error = np.zeros(numNx)
    for j in range(numNx):
        file = files[i][j]
        f = h5py.File(path+file, 'r')
        error[j] = f['error'][:][-1]
        f.close()
        print(error)
    plt.plot(Nx, error)
    
plt.legend(loc="upper right")
matplotlib.rcParams.update({'font.size': 16})
plt.grid()
plt.xscale('log')
plt.xlabel('Nx')
ylabel = r'$||\frac{\phi_i - \phi}{\phi}||_\infty$'
plt.ylabel(ylabel)




