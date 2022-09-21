# -*- coding: utf-8 -*-
import sys, os
sys.path.append(os.getcwd()+"/../")
from src.input_files.reeds_init import ReedsInit
from src.solvers.fixed_source.solvers import FixedSource
import matplotlib.pyplot as plt
from mpi4py import MPI
import numpy as np
import time
import cProfile

if __name__ == "__main__":
    comm = MPI.COMM_WORLD
    rank = comm.Get_rank()
    nproc = comm.Get_size()
    procname = MPI.Get_processor_name()
    
    N = 2**12
    Nx = 160
    G = 1
    generator = "halton"
    solver = "LGMRES"
    data = ReedsInit(N=N, Nx=Nx, generator=generator)
    start = time.time()
    maxit = 25
    tol = 1e-12
    phi = FixedSource(data,solver=solver, maxit=maxit, tol=tol, save_data=False)

    stop = time.time()
    if (rank==0):
        print("time: ",stop-start)
        plt.plot(range(Nx),phi)
