# -*- coding: utf-8 -*-

import sys
sys.path.append("../")
from src.init_files.reeds_init import ReedsInit
from src.functions.source_iteration import SourceIteration

if __name__ == "__main__":
    # initialize problem data
    N = 2**10
    data = ReedsInit(N=N, generator="halton")

    SI = SourceIteration(data)
    SI.max_iter = 20
    # run source iteration
    SI.Run()
