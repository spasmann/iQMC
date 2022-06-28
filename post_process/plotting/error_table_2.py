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
from functions.functions import ReduceFlux, RelError, AbsError, PlotLine

import sys, os
sys.path.append(os.getcwd()+"/../../")
from src.init_files.reeds_solution import reeds_mcdc_sol, reeds_sol, reeds_julia_sol
from src.init_files.mg_init import MultiGroupInit


Nvals = np.array((2**10, 2**11, 2**12,2**13,2**14)) #Reeds
problem="reeds_data"
NxBase = 80
#Nvals = np.array((2**7, 2**8, 2**9, 2**10,2**11)) #multigroup
#problem = "12"
#NxBase=10


generator="sobol"
nproc = 4
HaltonNx80 = PlotLine(Nvals=Nvals, Nx=NxBase, generator=generator, problem=problem, nproc=nproc)
HaltonNx160 = PlotLine(Nvals=Nvals, Nx=NxBase*2, generator=generator, problem=problem, nproc=nproc)
HaltonNx320 = PlotLine(Nvals=Nvals, Nx=NxBase*4, generator=generator, problem=problem, nproc=nproc)
O1 = HaltonNx320[0]*(Nvals[0]/Nvals)


generator="random"
nproc=2
RandomNx80 = PlotLine(Nvals=Nvals, Nx=NxBase, generator=generator, problem=problem, nproc=nproc)
RandomNx160 = PlotLine(Nvals=Nvals, Nx=NxBase*2, generator=generator, problem=problem, nproc=nproc)
RandomNx320 = PlotLine(Nvals=Nvals, Nx=NxBase*4, generator=generator, problem=problem, nproc=nproc)
O2 = RandomNx80[0]*np.sqrt(Nvals[0]/Nvals)

plt.figure(dpi=300, figsize=(10,5))
plt.suptitle("Reeds Problem Absolute Error")
ylabel = r'$||\phi - \phi_t||_\infty$'
ylim = [1e-2, 1e1]

plt.subplot(121)
plt.plot(Nvals, O1, label=r'$N^{-1}$')
plt.plot(Nvals, HaltonNx80,'o-',label=r'$N_x=80$')
plt.plot(Nvals, HaltonNx160,'s-',label=r'$N_x=160$')
plt.plot(Nvals, HaltonNx320,'^-',label=r'$N_x=320$')
plt.yscale('log')
plt.ylabel(ylabel)
plt.ylim(ylim)
plt.xlabel('N')
plt.xscale('log')
plt.title('QMC')
plt.legend()
plt.tight_layout()

plt.subplot(122)
plt.plot(Nvals, O2, label=r'$N^{-1/2}$')
plt.plot(Nvals, RandomNx80,'o-',label=r'$N_x=80$')
plt.plot(Nvals, RandomNx160,'s-',label=r'$N_x=160$')
plt.plot(Nvals, RandomNx320,'^-',label=r'$N_x=320$')
plt.yscale('log')
plt.ylabel(ylabel)
plt.ylim(ylim)
plt.xlabel('N')
plt.xscale('log')
plt.title('MC')
plt.legend()
plt.tight_layout()








