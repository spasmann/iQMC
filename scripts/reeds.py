# -*- coding: utf-8 -*-

import sys
sys.path.append("../")
from src.init_files.reeds_init import ReedsInit
from src.functions.source_iteration import SourceIteration

from time import process_time

if __name__ == "__main__":
    # initialize problem data
    N = 2**10
    Nx = 1
    start = process_time()
    data = ReedsInit(N=N, Nx=Nx, generator="halton")
    stop = process_time()
    time1 = stop-start
    print("SI Elapsed Time: {}s".format(round(time1,ndigits=5)))
    
    N = 2**10
    Nx = 180*10
    start = process_time()
    data = ReedsInit(N=N, Nx=Nx, generator="halton")
    stop = process_time()
    time2 = stop-start
    print("SI Elapsed Time: {}s".format(round(time2,ndigits=5)))
    print("A {}x speed up".format(round(time1/time2)))
    
    #SI = SourceIteration(data)
    #SI.max_iter = 10
    #SI.Run()

    

    