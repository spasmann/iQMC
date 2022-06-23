# -*- coding: utf-8 -*-
import sys
sys.path.append("../")
from src.init_files.reeds_init import ReedsInit
from src.init_files.reeds_solution import reeds_julia_sol
from src.solvers.solvers import LGMRES
import matplotlib.pyplot as plt
from mpi4py import MPI
import numpy as np
import time

if __name__ == "__main__":
    comm = MPI.COMM_WORLD
    rank = comm.Get_rank()
    nproc = comm.Get_size()
    procname = MPI.Get_processor_name()
    
    N = 2**12
    Nx = 80
    G = 1
    generator = "halton"
    data = ReedsInit(N=N, Nx=Nx, generator=generator)
    start = time.time()
    phi = LGMRES(data)
    stop = time.time()
    if (rank == 0):
        print("Time: ",stop-start)
        julia = reeds_julia_sol()
        print("Julia-Python Diff: ")
        print(abs(julia - phi).max())
        

    
        
        