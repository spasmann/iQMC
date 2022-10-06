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
from src.solvers.eigenvalue.davidson import Davidson
import matplotlib.pyplot as plt
import time

if __name__ == "__main__":
    # initialize problem data
    Nx = 20
    N = 2**12
    solver = "GMRES"
    generator = "halton"
    
    data = PUa_1_0_SL_init(N=N, Nx=Nx, generator=generator)
    data.save_data = False
    start = time.time()
    phi, khist = PowerIteration(data,
                        solver=solver,
                        max_outter_itt=25, 
                        max_inner_itt=25, 
                        outter_tol=1e-6,
                        inner_tol=1e-6)
    stop = time.time()
    
    data = PUa_1_0_SL_init(N=N, Nx=Nx, generator=generator)
    data.save_data = False
    start1 = time.time()
    phi2, keff, itt = Davidson(data, k0=1.0, l=1, m=4, numSweeps=8, tol=1e-6, maxit=25)
    stop1 = time.time()
    
    plt.figure(dpi=300)
    plt.plot(range(Nx),phi/phi.sum(), label='Power Iteration')
    plt.plot(range(Nx),phi2/phi2.sum(), '--', label='Davidson')
    plt.legend()
    plt.ylabel('Scalar Flux')
    plt.xlabel('Spatial Cells (x)')
    print("PI runtime: ", stop-start)
    print("Davidson runtime: ", stop1-start1)
    
    