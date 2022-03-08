# -*- coding: utf-8 -*-

import sys
sys.path.append("../")
from src.init_files.mg_init import MultiGroupInit
from src.functions.source_iteration import SourceIteration

if __name__ == "__main__":
    # initialize problem data
    N = 2**10
    #data1 = MultiGroupInit(N=N, generator="random")
    data2 = MultiGroupInit(N=N, generator="sobol")
    
    # initialize source iteration
    #SI = SourceIteration(data1)
    #SI.max_iter = 20
    # run source iteration
    #SI.Run()

    SI = SourceIteration(data2)
    SI.max_iter = 20
    # run source iteration
    SI.Run()


