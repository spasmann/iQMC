#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Apr  7 13:31:15 2022

@author: sampasmann
"""

import sys
sys.path.append("../")
from src.input_files.U_2_0_SL_init import U_2_0_SL_init
from src.solvers.eigenvalue.solvers import PowerIteration
from src.solvers.eigenvalue.davidson import Davidson
import matplotlib.pyplot as plt

if __name__ == "__main__":
    # initialize problem data
    Nx = 10
    N = 2**10
    maxit = 10
    generator = "halton"
    solver = "LGMRES"
    data = U_2_0_SL_init(N=N, Nx=Nx, generator=generator)
    data.save_data = False
    phi = PowerIteration(data,
                        solver=solver,
                        max_outter_itt=maxit, 
                        max_inner_itt=maxit, 
                        outter_tol=1e-6,
                        inner_tol=1e-6)
    #phi = Davidson(data, tol=1e-5, maxit=10)
  
    plt.plot(range(Nx),phi[0])
