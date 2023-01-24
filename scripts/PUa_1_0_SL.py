#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon Apr  4 17:52:16 2022

@author: sampasmann
"""

import sys
sys.path.append("../")
from src.input_files.PUa_1_0_SL_init import PUa_1_0_SL_init
from src.solvers.eigenvalue.solvers import PowerIteration, Davidson
import matplotlib.pyplot as plt
import time
import numpy as np

if __name__ == "__main__":
    # initialize problem data
    Nx          = 12
    N           = 2**10
    tol         = 1e-3
    maxit       = 25
    solver      = "LGMRES"
    generator   = "sobol"
    data        = PUa_1_0_SL_init(N=N, Nx=Nx, generator=generator, source_tilt=False)
    start = time.time()
    # phi, khist, itt = PowerIteration(data,
    #                     solver=solver,
    #                     max_outter_itt=maxit, 
    #                     max_inner_itt=maxit, 
    #                     outter_tol=tol,
    #                     inner_tol=tol)
    phi, khist, itt = Davidson(data, k0=1.0, l=1, m=None, numSweeps=4, tol=tol, maxit=30)
    stop = time.time()
    print("PI took: ", stop-start)
    plt.plot(range(Nx),data.tallies.phi_avg)
    
    
    # =============================================================================
    # plot source
    # =============================================================================
    source_tilt = data.source_tilt

    plt.figure(figsize=(6,4),dpi=300)
    plt.title('Source')
    q    = data.tallies.q
    mesh = data.mesh.edges
    mid  = data.mesh.midpoints
    x    = np.linspace(data.LB, data.RB, num=1000)
    n    = len(x)
    conditions = [(mesh[i] <= x) & (x <= mesh[i+1]) for i in range(Nx)]
    y1 = np.piecewise(x, conditions, q)
    plt.plot(x,y1, label=r'$a_j$')
    if (source_tilt):
        qdot        = data.tallies.qdot
        y2          = np.zeros_like(x)
        for i in range(n):
            zone = data.mesh.GetZone([x[i],0,0], [0,0,0])
            x_mid = data.mesh.midpoints[zone]
            y2[i] = q[zone] + qdot[zone]*(x[i] - x_mid)
        plt.plot(x,y2,label=r'$a_j + b_j(x)$')
    for i in range(len(mesh)):
        plt.axvline(mesh[i],linestyle='-',color='black')
    plt.legend()
    plt.tight_layout()
    
    #plt.plot(range(len(phi_hist)), phi_hist)
    #plt.yscale('log')
