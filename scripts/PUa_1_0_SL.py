#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon Apr  4 17:52:16 2022

@author: sampasmann
"""

import sys
sys.path.append("../")
from src.input_files.PUa_1_0_SL_init import PUa_1_0_SL_init
from src.solvers.eigenvalue.solvers import PowerIteration
import matplotlib.pyplot as plt

if __name__ == "__main__":
        # initialize problem data
    Nx = 20
    N = 2**12
    solver = "LGMRES"
    
    generator = "halton"
    data = PUa_1_0_SL_init(N=N, Nx=Nx, generator=generator)
    data.save_data = False
    phi, khist_halton = PowerIteration(data,
                        solver=solver,
                        max_outter_itt=25, 
                        max_inner_itt=25, 
                        outter_tol=1e-6,
                        inner_tol=1e-6)
    
    generator = "random"
    data = PUa_1_0_SL_init(N=N, Nx=Nx, generator=generator)
    data.save_data = False
    phi, khist_rand = PowerIteration(data,
                        solver=solver,
                        max_outter_itt=25, 
                        max_inner_itt=25, 
                        outter_tol=1e-9,
                        inner_tol=1e-9)
    
    plt.figure(dpi=300)
    plt.plot(range(khist_halton.size), abs(khist_halton-1), label="QMC")
    plt.plot(range(khist_rand.size), abs(khist_rand-1), label="MC")
    plt.yscale('log')
    plt.ylabel('K-Effective Abs Error')
    plt.xlabel('Iteration')
    plt.legend()
