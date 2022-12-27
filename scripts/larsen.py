# -*- coding: utf-8 -*-
import sys, os
sys.path.append(os.getcwd()+"/../")
from src.input_files.larsen_init import LarsenInit
from src.solvers.fixed_source.solvers import FixedSource
import matplotlib.pyplot as plt
import time
import numpy as np


if __name__ == "__main__":
    
    N           = 2**11
    Nx          = 10
    generator   = "sobol"
    solver      = "GMRES"
    data1       = LarsenInit(N=N, Nx=Nx, generator=generator, source_tilt=True)
    start       = time.time()
    maxit       = 1
    tol         = 1e-4
    phi1         = FixedSource(data1,solver=solver, maxit=maxit, tol=tol, 
                              save_data=False)
    sol = data1.true_flux
    stop = time.time()
    print("time: ",stop-start)
    plt.plot(data1.mesh.midpoints, phi1[:,0], label='iQMC')
    plt.plot(data1.mesh.midpoints, sol, label='Sol')
    plt.legend()