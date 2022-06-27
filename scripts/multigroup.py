# -*- coding: utf-8 -*-

import sys, os
sys.path.append(os.getcwd()+"/../")
from src.init_files.mg_init import MultiGroupInit
from src.solvers.solvers import LGMRES
import matplotlib.pyplot as plt
from mpi4py import MPI
import time


if __name__ == "__main__":
    comm = MPI.COMM_WORLD
    rank = comm.Get_rank()
    nproc = comm.Get_size()
    procname = MPI.Get_processor_name()
    
    N = 2**12
    Nx = 80
    G = 12
    generator = "halton"
    data = MultiGroupInit(numGroups=G, N=N, Nx=Nx, generator=generator)
    start = time.time()
    phi = LGMRES(data)
    stop = time.time()
    if (rank == 0):
        print("Time: ",stop-start)
        sol = data.true_flux
        print("Sol-QMC Diff: ")
        print(abs(sol - phi).max())
        