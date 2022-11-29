# -*- coding: utf-8 -*-
import sys, os
sys.path.append(os.getcwd()+"/../")
from src.input_files.garcia_init import GarciaInit
from src.solvers.fixed_source.solvers import FixedSource
import matplotlib.pyplot as plt
from mpi4py import MPI
import time


if __name__ == "__main__":
    comm = MPI.COMM_WORLD
    rank = comm.Get_rank()
    nproc = comm.Get_size()
    procname = MPI.Get_processor_name()
    # initialize problem data
    N           = 2**8
    Nx          = 60
    G           = 1
    generator   = "random"
    solver      = "Picard"
    data        = GarciaInit(N=N, Nx=Nx, generator=generator, source_tilt=True)
    start       = time.time()
    maxit       = 10
    tol         = 1e-4
    phi         = FixedSource(data,solver=solver, maxit=maxit, tol=tol, 
                              save_data=False)

    stop = time.time()
    if (rank==0):
        print("time: ",stop-start)
        plt.plot(range(Nx),phi)
