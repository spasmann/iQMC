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

if __name__ == "__main__":
    # initialize problem data
    Nx = 140
    N = 2**14
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
    
    if (rank==0):
        plt.plot(range(Nx),phi)