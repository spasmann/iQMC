# -*- coding: utf-8 -*-
"""
Created on Thu Jan 12 16:22:00 2023

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
    N_list      = 100*np.array((1,2,3,4,5,6,7,8,9,10,11,12,13,14,15,16,17,18,19,20))
    # Nx_list     = [2, 4, 8, 16, 32]
    Nx          = 4
    maxit       = 50
    tol         = 1e-9
    solver      = "LGMRES"
    generator   = "halton"
    
    err1        = np.zeros(N_list.size)
    err2        = np.zeros(N_list.size)
    
    err3        = np.zeros(N_list.size)
    err4        = np.zeros(N_list.size)
    
    for i in range(N_list.size):
        N           = N_list[i]
        # Nx          = Nx_list[i]
        
        data1       = LarsenInit(N=N, Nx=Nx, generator="halton",source_tilt=False)
        data2       = LarsenInit(N=N, Nx=Nx, generator="random",source_tilt=False)
        
        data3       = LarsenInit(N=N, Nx=Nx, generator="halton",source_tilt=True)
        data4       = LarsenInit(N=N, Nx=Nx, generator="random",source_tilt=True)
        
        phi1        = FixedSource(data1,solver=solver, maxit=maxit, tol=tol)
        phi2        = FixedSource(data2,solver=solver, maxit=maxit, tol=tol)
        
        phi3        = FixedSource(data3,solver=solver, maxit=maxit, tol=tol)
        phi4        = FixedSource(data4,solver=solver, maxit=maxit, tol=tol)
        
        err1[i]     = np.linalg.norm((phi1-data1.true_flux)/data1.true_flux.sum())
        err2[i]     = np.linalg.norm((phi2-data2.true_flux)/data2.true_flux.sum())
        
        err3[i]     = np.linalg.norm((phi3-data3.true_flux)/data3.true_flux.sum())
        err4[i]     = np.linalg.norm((phi4-data4.true_flux)/data4.true_flux.sum())

# =============================================================================
# 
# =============================================================================
        
        
plt.figure(dpi=300,figsize=(8,5))
plt.title('Cell Avg Scalar Flux Relative Error')

plt.plot(N_list, err1[0]*N_list[0]/N_list, label='$N^{-1}$')
plt.plot(N_list, err2[0]*np.sqrt(N_list[0]/N_list), label=r'$N^{-1/2}$')

plt.plot(N_list, err1, 'b-o', label='QMC')#label=r'$a_j$')
plt.plot(N_list, err2, 'g-^', label='MC')#label=r'$a_j + b_j(x)$')

plt.plot(N_list, err3, 'b--o', label='QMC (ST)')#label=r'$a_j$')
plt.plot(N_list, err4, 'g--^', label='MC (ST)')#label=r'$a_j + b_j(x)$')


plt.yscale('log')
plt.xscale('log')
plt.legend()


# =============================================================================
# subplots grouped by MC/QMC
# =============================================================================

plt.subplots(1,2, sharex=False, sharey=True, figsize=(9,4), dpi=300)
plt.suptitle('Larsen et al. Infinite Linear Source Problem')
count = 10

plt.subplot(121)
plt.title('QMC')
plt.plot(N_list[:count], err1[0]*N_list[0]/N_list[:count], 'k-',label='$N^{-1}$')
plt.plot(N_list[:count], err3[0]*N_list[0]/N_list[:count], 'k-')
plt.plot(N_list[:count], err1[:count], 'b-o', label='Constant Source')
plt.plot(N_list[:count], err3[:count], 'g--^', label='Linear Source')
plt.yscale('log')
plt.xscale('log')
plt.ylabel(r'$\phi$ Relative Error')
plt.xlabel('N')
plt.legend()

plt.subplot(122)
plt.title('MC')
plt.plot(N_list[:count], err2[0]*np.sqrt(N_list[0]/N_list[:count]), 'k-',label=r'$N^{-1/2}$')
plt.plot(N_list[:count], err2[:count], 'b-o', label='Constant Source')
plt.plot(N_list[:count], err4[:count], 'g--^', label='Linear Source')
plt.yscale('log')
plt.xscale('log')
plt.xlabel('N')
plt.legend()

# =============================================================================
# subplots grouped by Constant/Linear Source
# =============================================================================


plt.subplots(1,2, sharex=False, sharey=True, figsize=(9,4), dpi=300)
plt.suptitle('Larsen et al. Infinite Linear Source Problem')
count = 10

plt.subplot(121)
plt.title('Constant Source')
plt.plot(N_list[:count], err1[0]*N_list[0]/N_list[:count], 'k-',label='$N^{-1}$')
plt.plot(N_list[:count], err2[0]*np.sqrt(N_list[0]/N_list[:count]), 'k--',label='$N^{-1/2}$')
plt.plot(N_list[:count], err1[:count], 'b-o', label='QMC')
plt.plot(N_list[:count], err2[:count], 'g--^', label='MC')
plt.yscale('log')
plt.xscale('log')
plt.ylabel(r'$\phi$ Relative Error')
plt.xlabel('N')
plt.legend()

plt.subplot(122)
plt.title('Linear Source')
plt.plot(N_list[:count], err3[0]*N_list[0]/N_list[:count], 'k-',label='$N^{-1}$')
plt.plot(N_list[:count], err4[0]*np.sqrt(N_list[0]/N_list[:count]), 'k--',label='$N^{-1/2}$')
plt.plot(N_list[:count], err3[:count], 'b-o', label='QMC')
plt.plot(N_list[:count], err4[:count], 'g--^', label='MC')
plt.yscale('log')
plt.xscale('log')
plt.xlabel('N')
plt.legend()
