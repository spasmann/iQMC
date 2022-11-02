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
from src.solvers.eigenvalue.davidson import Davidson
import matplotlib.pyplot as plt

if __name__ == "__main__":
    # initialize problem data
    NList   = [2**7, 2**8, 2**9, 2**10, 2**11, 2**12, 2**13]
    Nx      = 20
    tol     = 1e-6
    maxit   = 25
    preC    = 6
    solver  = "LGMRES"
    QMC     = "halton"
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
        QMCphi, QMCkeff, QMCitt = Davidson(data, k0=1.0, l=1, m=maxit, numSweeps=preC, tol=tol, maxit=maxit)
        k_QMC_list.append(QMCkeff)
        
        data = PUa_1_0_SL_init(N=N, Nx=Nx, generator=MC)
        data.save_data = False
        MCphi, MCkeff, MCitt = Davidson(data, k0=1.0, l=1, m=maxit, numSweeps=preC, tol=tol, maxit=maxit)
        k_MC_list.append(MCkeff)
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
plt.plot(range(Nx), QMCphi, label="QMC")
plt.plot(range(Nx), MCphi, label="MC")
plt.legend()
plt.ylabel('Scalar Flux')
plt.xlabel('x')
plt.title('PUa-1-0-SL Scalar Flux')
