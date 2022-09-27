#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Apr  7 13:31:15 2022

@author: sampasmann
"""

import sys
sys.path.append("../")
from src.input_files.URRa_2_0_SL_init import URRa_2_0_SL_init
from src.solvers.eigenvalue.solvers import PowerIteration

if __name__ == "__main__":
    # initialize problem data
    Nx = 40
    N = 2**12
    generator = "halton"
    solver = "LGMRES"
    data = URRa_2_0_SL_init(N=N, Nx=Nx, generator=generator)
    data.save_data = False
    phi = PowerIteration(data,
                        solver=solver,
                        max_outter_itt=10, 
                        max_inner_itt=10, 
                        outter_tol=1e-6,
                        inner_tol=1e-6)
    
    if (rank==0):
        plt.plot(range(Nx),phi[0])
