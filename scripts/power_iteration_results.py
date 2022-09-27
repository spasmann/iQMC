#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon Sep 26 14:25:02 2022

@author: sampasmann
"""

import sys
sys.path.append("../")
from src.input_files.PUa_1_0_SL_init import PUa_1_0_SL_init
from src.solvers.eigenvalue.solvers import PowerIteration
import matplotlib.pyplot as plt

if __name__ == "__main__":
    # initialize problem data
    NList   = [2**10, 2**11, 2**12, 2**13, 2**14, 2**15, 2**16]
    Nx      = 40
    solver  = "LGMRES"
    QMC     = "halton"
    MC      = "random"
    k_QMC_list = []
    k_MC_list = []
    
    for i in range(len(NList)):
        N = NList[i]
        data = PUa_1_0_SL_init(N=N, Nx=Nx, generator=QMC)
        data.save_data = False
        phi_QMC, k_hist_QMC = PowerIteration(data,
                    solver=solver,
                    max_outter_itt=25, 
                    max_inner_itt=25, 
                    outter_tol=1e-9,
                    inner_tol=1e-9)
        k_QMC_list.append(k_hist_QMC[-1])
        
        data = PUa_1_0_SL_init(N=N, Nx=Nx, generator=MC)
        data.save_data = False
        phi_MC, k_hist_MC = PowerIteration(data,
                    solver=solver,
                    max_outter_itt=25, 
                    max_inner_itt=25, 
                    outter_tol=1e-9,
                    inner_tol=1e-9)
        k_MC_list.append(k_hist_MC[-1])


QMC = np.array(k_QMC_list)
MC = np.array(k_MC_list)
plt.figure(dpi=300)
plt.plot(NList, abs(1.0-QMC), label="QMC")
plt.plot(NList, abs(1.0-MC), label="MC")
plt.yscale('log')
plt.ylabel('Error')
plt.xscale('log')
plt.xlabel('N')
plt.legend()
plt.title('K-Effective Error')

plt.figure(dpi=300)
plt.plot(range(Nx), phi_QMC, label="QMC")
plt.plot(range(Nx), phi_MC, label="MC")
plt.legend()
plt.ylabel('Scalar Flux')
plt.xlabel('x')
plt.title('PUa-1-0-SL Scalar Flux')
