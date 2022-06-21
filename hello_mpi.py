# -*- coding: utf-8 -*-

from mpi4py import MPI

comm = MPI.COMM_WORLD
rank = comm.Get_rank()
size = comm.Get_size()
procname = MPI.Get_processor_name()

print("Hello world from processor {}, rank {:d} out of {:d} processors".format(procname, rank, size))