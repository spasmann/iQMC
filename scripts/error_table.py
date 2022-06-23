#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Jun 23 11:48:20 2022

@author: sampasmann
"""

import sys, os
sys.path.append(os.getcwd()+"/../")
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
    
    Nvals = np.array((2**10,2**11,2**12,2**13,2**14,2**15,2**16,2**17,2**18,2**19,2**20))
    Nxvals = 80*np.array((1,2,4))    
    generator = "halton"
    
    for Nx in Nxvals:
        for N in Nvals:
            data = ReedsInit(N=N,Nx=Nx,generator=generator)
            start = time.time()
            phi = LGMRES(data)
            stop = time.time()
            if (rank == 0):
                print("Nx: ", Nx, "N: ", N, "Time: ",stop-start)
            
            
            