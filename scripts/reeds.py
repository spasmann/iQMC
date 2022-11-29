# -*- coding: utf-8 -*-
import sys, os
sys.path.append(os.getcwd()+"/../")
from src.input_files.reeds_init import ReedsInit
from src.solvers.fixed_source.solvers import FixedSource
import matplotlib.pyplot as plt
from mpi4py import MPI
import time
import numpy as np

if __name__ == "__main__":
    comm = MPI.COMM_WORLD
    rank = comm.Get_rank()
    nproc = comm.Get_size()
    procname = MPI.Get_processor_name()
    
    N           = 2**11
    Nx          = 64
    generator   = "sobol"
    solver      = "LGMRES"
    data        = ReedsInit(N=N, Nx=Nx, generator=generator, source_tilt=True, LB=-8.0, RB=8.0)
    start       = time.time()
    maxit       = 10
    tol         = 1e-4
    phi2         = FixedSource(data,solver=solver, maxit=maxit, tol=tol, 
                              save_data=False)

    stop = time.time()
    if (rank==0):
        print("time: ",stop-start)
        plt.plot(data.mesh.midpoints, phi2)
        #plt.plot(data.mesh.midpoints, data.tallies.phi_avg[:,0] + data.tallies.dphi_s[:,0]*data.mesh.midpoints)
        #plt.plot(data.mesh.midpoints, np.ones(Nx), 'o')
