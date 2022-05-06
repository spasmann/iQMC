# -*- coding: utf-8 -*-


from src.init_files.first_init import FirstInit
from src.functions.source_iteration import SourceIteration

if __name__ == "__main__":
    # initialize problem data
    N = 2**10
    data1 = FirstInit(N=N, generator="random")
    data2 = FirstInit(N=N, generator="sobol")
    # initialize source iteration
    SI = SourceIteration(data1)
    SI.max_iter = 20
    # run source iteration
    SI.Run()
    
    SI = SourceIteration(data2)
    SI.max_iter = 20
    # run source iteration
    SI.Run()

