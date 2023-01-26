# -*- coding: utf-8 -*-
"""
Created on Wed Jan 25 16:03:23 2023

@author: Sam
"""

import sys, os
sys.path.append(os.getcwd()+"/../")
from src.input_files.garcia_init import GarciaInit
from src.solvers.fixed_source.solvers import FixedSource
from post_process.functions.functions import SN_Sweep, garcia_angular_flux_sol, garcia_angle_bins
import matplotlib.pyplot as plt
import time
import numpy as np

if __name__ == "__main__":
    N_list      = 100*np.array((1,2,3,4,5,6,7,8,9,10))
    # Nx_list     = [2, 4, 8, 16, 32]
    Nx          = 24
    maxit       = 50
    tol         = 1e-6
    s           = 1
    solver      = "LGMRES"
    generator   = "halton"
    sol         = garcia_angular_flux_sol(s)
    angles      = garcia_angle_bins()
    Na2         = angles.size
    
    err1        = np.zeros(N_list.size)
    err2        = np.zeros(N_list.size)
    err3        = np.zeros(N_list.size)
    err4        = np.zeros(N_list.size)
    
    
    for i in range(N_list.size):
        N           = N_list[i]
        # generate qmc data
        data1       = GarciaInit(N=N, Nx=Nx, generator="halton",source_tilt=False)
        data2       = GarciaInit(N=N, Nx=Nx, generator="random",source_tilt=False)
        data3       = GarciaInit(N=N, Nx=Nx, generator="halton",source_tilt=True)
        data4       = GarciaInit(N=N, Nx=Nx, generator="random",source_tilt=True)
        # execute iQMC
        phi1        = FixedSource(data1,solver=solver, maxit=maxit, tol=tol)
        phi2        = FixedSource(data2,solver=solver, maxit=maxit, tol=tol)
        phi3        = FixedSource(data3,solver=solver, maxit=maxit, tol=tol)
        phi4        = FixedSource(data4,solver=solver, maxit=maxit, tol=tol)
        # perform SN sweep to retrieve angular flux
        psi1        = SN_Sweep(angles, data1)
        psi2        = SN_Sweep(angles, data2)
        psi3        = SN_Sweep(angles, data3)
        psi4        = SN_Sweep(angles, data4)
        # relative error
        err1[i]     = np.linalg.norm((psi1-sol)/sol.sum())
        err2[i]     = np.linalg.norm((psi2-sol)/sol.sum())
        err3[i]     = np.linalg.norm((psi3-sol)/sol.sum())
        err4[i]     = np.linalg.norm((psi4-sol)/sol.sum())

# =============================================================================
# Error plot
# =============================================================================
        
        
# plt.figure(dpi=300,figsize=(8,5))
# plt.title('Cell Avg Scalar Flux Relative Error')

# plt.plot(N_list, err1[0]*N_list[0]/N_list, label='$N^{-1}$')
# plt.plot(N_list, err2[0]*np.sqrt(N_list[0]/N_list), label=r'$N^{-1/2}$')

# plt.plot(N_list, err1, 'b-o', label='QMC')
# plt.plot(N_list, err2, 'g-^', label='MC')

# plt.plot(N_list, err3, 'b--o', label='QMC (ST)')
# plt.plot(N_list, err4, 'g--^', label='MC (ST)')


# plt.yscale('log')
# plt.xscale('log')
# plt.legend()


# =============================================================================
# subplots grouped by Constant/Linear Source
# =============================================================================


plt.subplots(sharex=False, sharey=True, nrows=1,ncols=2,dpi=300,figsize=(8,4))
plt.suptitle('Garcia et al. Angular Flux Relative Error')
count = -1

plt.subplot(121)
plt.title('Constant Source')
plt.plot(N_list[:count], err1[0]*N_list[0]/N_list[:count], 'k-',label='$N^{-1}$')
plt.plot(N_list[:count], err2[0]*np.sqrt(N_list[0]/N_list[:count]), 'k--',label='$N^{-1/2}$')
plt.plot(N_list[:count], err1[:count], 'b-o', label='QMC')
plt.plot(N_list[:count], err2[:count], 'g--^', label='MC')
plt.yscale('log')
plt.xscale('log')
plt.ylabel(r'$\psi$ Relative Error')
plt.xlabel('N')

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

plt.tight_layout()