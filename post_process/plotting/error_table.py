#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue May 10 11:32:37 2022

@author: sampasmann
"""

import matplotlib.pyplot as plt
import matplotlib
import numpy as np
import h5py
import itertools
marker = itertools.cycle(('*', 's', '^', 'o', 'P','x','D')) 

path = "../saved_data/"
generator = "halton"
problem = "reeds_data"
N = [2**10, 2**11, 2**12, 2**13, 2**14, 2**15]
#Nx = [16*4, 16*8, 16*16, 16*32, 16*64, 16*128] # increasing spatial cells too fast
Nx = [16*4, 16*6, 16*8, 16*10, 16*12, 16*14, 16*16, 16*32, 16*128]  # second attempt 
Nx = [16*4, 16*6, 16*8, 16*10, 16*14, 16*16, 16*32] 
#Nx = [16*1, 16*2, 16*3, 16*4, 16*5, 16*6] # third attempt

N = np.array(N)
Nx = np.array(Nx)

plt.figure(dpi=300)
for i in range(len(Nx)):
    line = np.zeros(len(N))
    for j in range(len(N)):
        try:
            file = path+problem+"-"+generator+"-"+str(N[j])+"-"+str(Nx[i])
            f = h5py.File(file, 'r')    
            line[j] = f['error'][:][-1]
            f.close()
        except:
            print("File DNE: "+file)
        
    plt.plot(N/Nx[i],line,'o',label=str(Nx[i])) 

plt.yscale('log')
plt.xscale('log')
plt.xlabel("N/Nx")
ylabel = r'$||\phi - \phi_t||_\infty$'
plt.ylabel(ylabel)
plt.legend()
plt.title("Reeds Problem Absolute Error")
