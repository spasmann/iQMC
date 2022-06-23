# -*- coding: utf-8 -*-

from mpi4py import MPI
import numpy as np
from src.init_files.reeds_init import ReedsInit
from src.functions.source import GetSource
from src.functions.sweep import Sweep
from src.functions.tallies import Tallies


comm = MPI.COMM_WORLD
rank = comm.Get_rank()
nproc = comm.Get_size()
procname = MPI.Get_processor_name()

print("Hello world from processor {}, rank {:d} out of {:d} processors".format(procname, rank, nproc))
N = 2**10
Nx = 16
G = 1
    
vector = np.ones((2,2),dtype=float)*rank

vector = comm.allreduce(vector,op=MPI.SUM)
    
if (rank==0):
    print(vector)