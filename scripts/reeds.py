# -*- coding: utf-8 -*-
import sys, os
sys.path.append(os.getcwd()+"/../")
from src.init_files.reeds_init import ReedsInit
#from src.init_files.reeds_solution import reeds_julia_sol, reeds_mcdc_sol, reeds_sol
from src.solvers.solvers import LGMRES, Picard
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
    
    N = 2**14
    Nx = 80
    G = 1
    generator = "sobol"
    data = ReedsInit(N=N, Nx=Nx, generator=generator)
    start = time.time()
    maxit = 20
    phi = Picard(data,maxit=maxit, save_data=False,report_progress=True)
    stop = time.time()
    if (rank==0):
        print("time: ",stop-start)
        plt.plot(range(Nx),phi)

        

    
        
        
