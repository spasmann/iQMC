#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Jun 23 11:48:20 2022

@author: sampasmann
"""

import sys, os
sys.path.append(os.getcwd()+"/../")
from src.input_files.reeds_init import ReedsInit
from src.input_files.reeds_solution import reeds_sol
from src.solvers.fixed_source.solvers import FixedSource
import matplotlib.pyplot as plt
from mpi4py import MPI
import numpy as np
import time

if __name__ == "__main__":
    comm = MPI.COMM_WORLD
    rank = comm.Get_rank()
    nproc = comm.Get_size()
    procname = MPI.Get_processor_name()
    
    Nvals = np.array((2**7, 2**8, 2**9))
    #Nvals = np.array((2**10, 2**11, 2**12,2**13,2**14)) #Reeds
    Nx = 64
    generator = "sobol"
    solver = "LGMRES"
    maxit  = 15
    tol    = 1-6
    sol     = reeds_sol(Nx=Nx, LB=-8.0, RB=8.0)
    
    for N in Nvals:
        # Spatial Averaging
        data    = ReedsInit(N=N,Nx=Nx,generator=generator,source_tilt=False)
        start   = time.time()
        phi     = FixedSource(data,solver=solver, maxit=maxit, tol=tol,report_progress=False)
        stop    = time.time()
        if (rank == 0):
            print("Nx:", Nx, " N:", N, " Time:",stop-start)
            
            
            