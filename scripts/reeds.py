# -*- coding: utf-8 -*-
import sys, os
sys.path.append(os.getcwd()+"/../")
from src.input_files.reeds_init import ReedsInit
from src.solvers.fixed_source.solvers import FixedSource
import matplotlib.pyplot as plt
from mpi4py import MPI
import time
import numpy as np

from src.input_files.reeds_solution import reeds_mcdc_sol, reeds_sol

if __name__ == "__main__":
    comm = MPI.COMM_WORLD
    rank = comm.Get_rank()
    nproc = comm.Get_size()
    procname = MPI.Get_processor_name()
    
    N           = 2**11
    Nx          = 32
    generator   = "sobol"
    solver      = "LGMRES"
    data1       = ReedsInit(N=N, Nx=Nx, generator=generator, source_tilt=False)
    data2       = ReedsInit(N=N, Nx=Nx, generator=generator, source_tilt=True)
    start       = time.time()
    maxit       = 10
    tol         = 1e-4
    phi1         = FixedSource(data1,solver=solver, maxit=maxit, tol=tol, 
                              save_data=False)
    phi2         = FixedSource(data2,solver=solver, maxit=maxit, tol=tol, 
                              save_data=False)
    
    sol = reeds_mcdc_sol()
    sol = reeds_sol(Nx=32)
    
    err1 = np.linalg.norm(sol-phi1)
    err2 = np.linalg.norm(sol-phi2)
    print("Constant Source Error: ",err1)
    print("Linear Source Error: ",err2)

    stop = time.time()
    if (rank==0):
        print("time: ",stop-start)
        plt.plot(data1.mesh.midpoints, phi1)
        plt.plot(data1.mesh.midpoints, sol)