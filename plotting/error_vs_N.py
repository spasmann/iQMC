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
import itertools
marker = itertools.cycle(('*', 's', '^', 'o', 'P','x','D')) 

path = "../saved_data/"
generator = "halton"
problem = "reeds_data"
N = [2**9, 2**10, 2**11, 2**12, 2**13]
#Nx = [16*4, 16*8, 16*16, 16*32] # increasing spatial cells too fast
#Nx = [16*4, 16*6, 16*8, 16*10, 16*14]  # second attempt 
#Nx = [16*1, 16*2, 16*3, 16*4, 16*5, 16*6] # third attempt
Nx = [16*2, 16*4, 16*8, 16*16, 16*32] 
diagonal = True

plt.figure(dpi=300)
if (diagonal):
    diag = np.zeros(len(N))
for i in range(len(Nx)):
    line = np.zeros(len(N))
    for j in range(len(N)):
        try:
            file = path+problem+"-"+generator+"-"+str(N[j])+"-"+str(Nx[i])
            f = h5py.File(file, 'r')    
            line[j] = f['error'][:][-1]
            if (diagonal):
                if (i == j):
                   diag[j] = f['error'][:][-1]
            f.close()
    # theoretical convergence line
        except:
            print("File DNE: "+file)
    # which values of Nx to actually plot
    if (Nx[i] == 64) or (Nx[i] == 128) or (Nx[i] == 256):
        plt.plot(N,line,label="Nx="+str(Nx[i]),marker=next(marker))
    if (Nx[i] == 128):
        # 1/sqrt(N) convergence (MC)
        #plt.plot(N,line[0]*np.sqrt(N[0])/np.sqrt(N),'r-',label=r"$O(N^{-0.5})$")
        # 1/N convergence (QMC)
        plt.plot(N,line[0]*N[0]/N,'r-',label=r"$O(N^{-1})$")
    

if (diagonal):
    plt.plot(N,diag,"k--",label=r'$Nx(n)=16*2^{n}$')
plt.yscale('log')
plt.xscale('log')
plt.xlabel("N")
ylabel = r'$||\phi - \phi_t||_\infty$'
plt.ylabel(ylabel)
plt.legend()
plt.title("Reeds Problem Absolute Error")
