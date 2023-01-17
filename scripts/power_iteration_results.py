#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon Sep 26 14:25:02 2022

@author: sampasmann
"""

import sys
sys.path.append("../")
import numpy as np
from src.input_files.PUa_1_0_SL_init import PUa_1_0_SL_init
from src.solvers.eigenvalue.solvers import PowerIteration
import matplotlib.pyplot as plt

if __name__ == "__main__":
    # initialize problem data
    NList   = [2**7]
    Nx      = 20
    tol     = 1e-6
    maxit   = 25
    solver  = "LGMRES"
    QMC     = "sobol"
    MC      = "random"
    k_QMC_list = []
    k_MC_list = []
    print("    Keff Error    ")
    print("------------------")
    print("   MC        QMC   ")
    for i in range(len(NList)):
        N = NList[i]
        data = PUa_1_0_SL_init(N=N, Nx=Nx, generator=QMC)
        data.save_data = False
        phi_QMC, k_hist_QMC, itt_QMC = PowerIteration(data,
                                             solver         = solver,
                                             max_outter_itt = maxit, 
                                             max_inner_itt  = maxit, 
                                             outter_tol     = tol,
                                             inner_tol      = tol)
        k_QMC_list.append(k_hist_QMC[-1])
        
        data = PUa_1_0_SL_init(N=N, Nx=Nx, generator=MC)
        data.save_data = False
        phi_MC, k_hist_MC, itt_MC = PowerIteration(data,
                                           solver         = solver,
                                           max_outter_itt = maxit, 
                                           max_inner_itt  = maxit, 
                                           outter_tol     = tol,
                                           inner_tol      = tol)
        k_MC_list.append(k_hist_MC[-1])
        print(abs(1.0-k_MC_list[-1]), "   ", abs(1.0-k_QMC_list[-1]))


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
