# -*- coding: utf-8 -*-

import sys
sys.path.append("../")
from src.init_files.reeds_init import ReedsInit
from src.functions.source_iteration import SourceIteration

if __name__ == "__main__":
    # initialize problem data
    N = 2**11
    data = ReedsInit(N=N, generator="sobol")

    SI = SourceIteration(data)
    SI.max_iter = 10
    SI.Run()
    
    #data.moment_match = True
    #SI = SourceIteration(data)
    #SI.fname = "moment_match_halton"
    #SI.max_iter = 10
    # run source iteration
    #SI.Run()
