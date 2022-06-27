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

import sys, os
sys.path.append(os.getcwd()+"/../")
from functions.functions import ReduceFlux, RelError, AbsError
import sys, os
sys.path.append(os.getcwd()+"/../../")
from src.init_files.reeds_solution import reeds_mcdc_sol, reeds_sol, reeds_julia_sol

def PlotLine(Nvals=np.array((2**10)), Nx=80, generator="halton", nproc=64):
    sol = reeds_mcdc_sol()
    line = np.empty(len(Nvals))
    path = "../saved_data/"
    generator = "halton"
    problem = "reeds_data"
    
    for i in range(len(Nvals)):
        try:
            file = path+problem+"-"+generator+"-"+str(Nvals[i])+"-"+str(Nx)+"-"+str(nproc)
            f = h5py.File(file, 'r')    
            phi = f['phi_avg'][:]
            f.close()
            
            phi = ReduceFlux(phi, 80)
            error = RelError(phi, sol)
            line[i] = error
            
        except:
            print("File or Variable DNE: "+file)
            break

    return line


Nvals = np.array((2**10, 2**11, 2**12,2**13,2**14,2**15,2**16,2**17,2**18,2**19,2**20))

HaltonNx80 = PlotLine(Nvals=Nvals, Nx=80)
HaltonNx160 = PlotLine(Nvals=Nvals, Nx=160)
HaltonNx320 = PlotLine(Nvals=Nvals, Nx=320)
O = HaltonNx80[0]*(Nvals[0]/Nvals)

plt.figure(dpi=300, figsize=(10,5))
plt.suptitle("Reeds Problem Absolute Error")
ylabel = r'$||\phi - \phi_t||_\infty$'

plt.subplot(121)
plt.plot(Nvals, O, label=r'$N^{-1}$')
plt.plot(Nvals, HaltonNx80,'o-',label=r'$N_x=80$')
plt.plot(Nvals, HaltonNx160,'s-',label=r'$N_x=160$')
plt.plot(Nvals, HaltonNx320,'^-',label=r'$N_x=320$')
plt.yscale('log')
plt.ylabel(ylabel)
plt.xlabel('N')
plt.xscale('log')
plt.title('QMC')
plt.legend()
plt.tight_layout()

plt.subplot(122)
plt.yscale('log')
plt.ylabel(ylabel)
plt.xlabel('N')
plt.xscale('log')
plt.title('MC')
#plt.legend()
plt.tight_layout()








