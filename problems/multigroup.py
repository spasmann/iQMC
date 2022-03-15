# -*- coding: utf-8 -*-

import sys
sys.path.append("../")
from src.init_files.mg_init import MultiGroupInit
from src.functions.source_iteration import SourceIteration

if __name__ == "__main__":
    # initialize problem data
    N = 2**12
    data1 = MultiGroupInit(N=N, generator="random")
    data2 = MultiGroupInit(N=N, generator="halton")
    data3 = MultiGroupInit(N=N, generator="sobol")
    

    SI = SourceIteration(data1)
    SI.max_iter = 20
    SI.Run()
    SI = SourceIteration(data2)
    SI.max_iter = 20
    SI.Run()
    SI = SourceIteration(data3)
    SI.max_iter = 20
    SI.Run()


