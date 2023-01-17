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
    N_list      = np.array((2**10, 2**11, 2**12))
    # Nx_list     = [2, 4, 8, 16, 32]
    Nx          = 1
    maxit       = 50
    tol         = 1e-6
    solver      = "LGMRES"
    generator   = "random"
    
    err1        = np.zeros(N_list.size)
    err2        = np.zeros(N_list.size)
    
    for i in range(N_list.size):
        N           = N_list[i]
        # Nx          = Nx_list[i]
        
        data1       = LarsenInit(N=N, Nx=Nx, generator=generator,source_tilt=False)
        data2       = LarsenInit(N=N, Nx=Nx, generator=generator,source_tilt=True)
        
        phi1        = FixedSource(data1,solver=solver, maxit=maxit, tol=tol)
        phi2        = FixedSource(data2,solver=solver, maxit=maxit, tol=tol)
        
        err1[i]     = np.linalg.norm((phi1-data1.true_flux))
        err2[i]     = np.linalg.norm((phi2-data2.true_flux))
        
        plt.figure(1, dpi=300)
        plt.plot(data1.mesh.midpoints, data1.true_flux)
        plt.plot(data1.mesh.midpoints, phi1, '-o', label=r'$a_j$')
        plt.plot(data2.mesh.midpoints, phi2, '-o', label=r'$a_j + b_j(x)$')
        plt.legend()
        plt.show()
        
        
plt.figure(dpi=300)
plt.plot(N_list, err1[0]*N_list[0]/N_list, label='$N^{-1}$')
plt.plot(N_list, err1[0]*np.sqrt(N_list[0]/N_list), label=r'$N^{-1/2}$')
plt.plot(N_list, err1, '-o', label=r'$a_j$')
plt.plot(N_list, err2, '--^', label=r'$a_j + b_j(x)$')
plt.yscale('log')
plt.xscale('log')
plt.legend()
