#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon Apr  4 17:52:16 2022

@author: sampasmann
"""

import sys
sys.path.append("../")
from src.input_files.PUa_1_0_SP_init import PUa_1_0_SP_init
from src.solvers.eigenvalue.solvers import PowerIteration, Davidson
import matplotlib.pyplot as plt
import time

if __name__ == "__main__":
    # initialize problem data
    Nx = 20
    N = 2**11
    solver = "LGMRES"
    generator = "sobol"
    
    data = PUa_1_0_SP_init(N=N, Nx=Nx, generator=generator, source_tilt=True)
    start = time.time()
    #phi = Davidson(data, tol=1e-5, maxit=10)
    
    phi, khist, itt = PowerIteration(data,
                        solver=solver,
                        max_outter_itt=10, 
                        max_inner_itt=10, 
                        outter_tol=1e-4,
                        inner_tol=1e-4)
    stop = time.time()
    print("PI took: ", stop-start)
    plt.plot(data.mesh.midpoints, phi)
