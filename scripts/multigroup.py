# -*- coding: utf-8 -*-

import sys
sys.path.append("../")
from src.init_files.mg_init import MultiGroupInit
from src.functions.source_iteration import SourceIteration

import numpy as np

if __name__ == "__main__":
    # initialize problem data
    N = 2**10
    Nx = 5
    data1 = MultiGroupInit(N=N, Nx=Nx, generator="halton")
    SI = SourceIteration(data1)
    SI.max_iter = 50
    SI.Run()
    
