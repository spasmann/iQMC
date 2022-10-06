#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon Apr  4 18:29:16 2022

@author: sampasmann
"""

import sys
sys.path.append("../")
from src.input_files.PUa_H2O_1_0_SL_init import PUa_H2O_1_0_SL_init
from src.solvers.eigenvalue.solvers import PowerIteration
from src.solvers.eigenvalue.davidson import Davidson

if __name__ == "__main__":
    # initialize problem data
    Nx = 57
    N = 5000
    generator = "halton"
    solver = "LGMRES"
    data = PUa_H2O_1_0_SL_init(N=N, Nx=Nx, generator=generator)
    data.save_data = False
    phi = PowerIteration(data,
                        solver=solver,
                        max_outter_itt=25, 
                        max_inner_itt=25, 
                        outter_tol=1e-6,
                        inner_tol=1e-6)
    
    data = PUa_H2O_1_0_SL_init(N=N, Nx=Nx, generator=generator)
    data.save_data = False
    phi2, keff, itt = Davidson(data, k0=1.0, l=1, m=4, numSweeps=8, tol=1e-6, maxit=25)
    
    if (rank==0):
        plt.figure(dpi=300)
        plt.plot(range(Nx),phi[0], label='Power Iteration')
        plt.plot(range(Nx),phi2,label='Davidson')
        plt.legend()
        plt.ylabel('Scalar Flux')
        plt.xlabel('Spatial Cell (x)')