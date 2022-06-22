# -*- coding: utf-8 -*-
import sys
sys.path.append("../")
from src.init_files.reeds_init import ReedsInit
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
    
    N = 2**11
    Nx = 16*5
    G = 1
    generator = "halton"
    data = ReedsInit(N=N, Nx=Nx, generator=generator)
    start = time.time()
    phi = LGMRES(data)
    stop = time.time()
    if (rank == 0):
        print("Time: ",stop-start)
        sol = data.true_flux
        #sol = np.reshape(sol, (Nx*G,))
        #phi = np.reshape(phi, (Nx*G,))
        infnorm = np.linalg.norm(sol-phi)
        print("Error: ",infnorm)
        

    
        
        