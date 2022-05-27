# -*- coding: utf-8 -*-

import sys
sys.path.append("../")
from src.init_files.garcia_init import GarciaInit
from src.functions.source_iteration import SourceIteration

if __name__ == "__main__":
    # initialize problem data
    N = 2**9
    data = GarciaInit(N=N, generator="halton")
    # initialize source iteration
    SI = SourceIteration(data)
    SI.max_iter = 20
    # run source iteration
    SI.Run()


