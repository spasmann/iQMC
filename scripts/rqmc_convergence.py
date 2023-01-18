# -*- coding: utf-8 -*-
"""
Created on Wed Jan 18 09:12:52 2023

@author: Sam
"""

import sys, os
sys.path.append(os.getcwd()+"/../")
from src.input_files.larsen_init import LarsenInit
from src.solvers.fixed_source.solvers import FixedSource
import matplotlib.pyplot as plt
import time
import numpy as np

if __name__ == "__main__":
    N_list      = 100*np.array((1,2,3,4,5,6,7,8,9))
    J           = 40
    seeds       = np.array((1000000*np.random.random(J)),dtype=int)
    Nx          = 2
    maxit       = 50
    tol         = 1e-6
    solver      = "LGMRES"
    generator   = "halton"
    err1        = np.zeros((N_list.size, seeds.size))
    err2        = np.zeros((N_list.size, seeds.size))
    
    for i in range(N_list.size):
        N = N_list[i]
        for j in range(seeds.size):
            data1       = LarsenInit(N=N, Nx=Nx, generator="halton",
                                     source_tilt=False, RQMC=True, seed=seeds[j])
            phi1        = FixedSource(data1,solver=solver,maxit=maxit,tol=tol)
            err1[i,j]   = np.linalg.norm((phi1-data1.true_flux)/data1.true_flux.sum())
            
            data2       = LarsenInit(N=N, Nx=Nx, generator="halton",
                                     source_tilt=True, RQMC=True, seed=seeds[j])
            phi2        = FixedSource(data2,solver=solver,maxit=maxit,tol=tol)
            err2[i,j]   = np.linalg.norm((phi2-data2.true_flux)/data2.true_flux.sum())

# =============================================================================
# 
# =============================================================================

count1 = 2
count2 = N_list.size
        
plt.figure(dpi=300,figsize=(8,5))
plt.title('Cell Avg Scalar Flux Relative Error')

avg_err1 = err1.sum(axis=1)/seeds.size
avg_err2 = err2.sum(axis=1)/seeds.size

plt.plot(N_list[count1:count2], avg_err2[count1]*N_list[count1]/N_list[count1:count2], label='$N^{-1}$')
plt.plot(N_list[count1:count2], avg_err1[count1]*np.sqrt(N_list[count1]/N_list[count1:count2]), label='$N^{-1/2}$')

plt.plot(N_list[count1:count2], avg_err1[count1:count2], 'b-o', label='Linear')
plt.plot(N_list[count1:count2], avg_err2[count1:count2], 'g-^', label='Constant')

plt.yscale('log')
plt.xscale('log')
plt.legend()
            
# =============================================================================
# subplots grouped by Constant/Linear Source
# =============================================================================


# plt.subplots(1,2, sharex=False, sharey=True, figsize=(9,4), dpi=300)
# plt.suptitle('Larsen et al. Infinite Linear Source Problem')
# count = 10

# plt.subplot(121)
# plt.title('Constant Source')
# plt.plot(N_list[:count], err1[0]*N_list[0]/N_list[:count], 'k-',label='$N^{-1}$')
# plt.plot(N_list[:count], err2[0]*np.sqrt(N_list[0]/N_list[:count]), 'k--',label='$N^{-1/2}$')
# plt.plot(N_list[:count], err1[:count], 'b-o', label='QMC')
# plt.plot(N_list[:count], err2[:count], 'g--^', label='MC')
# plt.yscale('log')
# plt.xscale('log')
# plt.ylabel(r'$\phi$ Relative Error')
# plt.xlabel('N')
# plt.legend()

# plt.subplot(122)
# plt.title('Linear Source')
# plt.plot(N_list[:count], err3[0]*N_list[0]/N_list[:count], 'k-',label='$N^{-1}$')
# plt.plot(N_list[:count], err4[0]*np.sqrt(N_list[0]/N_list[:count]), 'k--',label='$N^{-1/2}$')
# plt.plot(N_list[:count], err3[:count], 'b-o', label='QMC')
# plt.plot(N_list[:count], err4[:count], 'g--^', label='MC')
# plt.yscale('log')
# plt.xscale('log')
# plt.xlabel('N')
# plt.legend()