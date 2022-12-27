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
    
    N           = 1
    Nx          = 32
    generator   = "halton"
    solver      = "Picard"
    data1       = ReedsInit(N=N, Nx=Nx, generator=generator, source_tilt=False)
    start       = time.time()
    maxit       = 1
    tol         = 1e-4
    phi1         = FixedSource(data1,solver=solver, maxit=maxit, tol=tol, 
                              save_data=False)
    sol = reeds_sol(Nx=Nx)
    stop = time.time()
    if (rank==0):
        print("time: ",stop-start)
        plt.plot(data1.mesh.midpoints, phi1[:,0], label='iQMC')
        # plt.plot(data1.mesh.midpoints, sol, label='Sol')
        plt.legend()