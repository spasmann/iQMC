# -*- coding: utf-8 -*-

import sys
sys.path.append("../")
from src.init_files.garcia_init import GarciaInit
from src.functions.source_iteration import SourceIteration
#from src.solvers.solvers import Picard

if __name__ == "__main__":
    # initialize problem data
    N = 2**11
    Nx = 402
    qmc_data = GarciaInit(N=N, Nx=Nx, generator="halton")
    # initialize source iteration
    #maxit = 30
    #data_out = Picard(qmc_data,maxit=maxit,report_progress=True)
    SI = SourceIteration(qmc_data)
    SI.max_iter = 30
    # run source iteration
    SI.Run()


