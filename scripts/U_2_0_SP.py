#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Apr  7 13:31:15 2022

@author: sampasmann
"""

import sys
sys.path.append("../")
from src.input_files.U_2_0_SP_init import U_2_0_SP_init
from src.solvers.eigenvalue.solvers import PowerIteration
from src.solvers.eigenvalue.davidson import Davidson
import matplotlib.pyplot as plt

if __name__ == "__main__":
    # initialize problem data
    Nx = 11
    N = 1000
    maxit = 10
    generator = "halton"
    solver = "LGMRES"
    data = U_2_0_SP_init(N=N, Nx=Nx, generator=generator)
    data.save_data = False
    phi = PowerIteration(data,
                        solver=solver,
                        max_outter_itt=maxit, 
                        max_inner_itt=maxit, 
                        outter_tol=1e-5,
                        inner_tol=1e-5)
    #phi = Davidson(data, tol=1e-5, maxit=10)
  
    plt.plot(range(Nx),phi[0])
