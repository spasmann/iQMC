# -*- coding: utf-8 -*-

import sys
sys.path.append("../")
from src.init_files.reeds_init import ReedsInit
from src.functions.source_iteration import SourceIteration

if __name__ == "__main__":
    N = 2**11
    Nx = 320
    data = ReedsInit(N=N, Nx=Nx, generator="latin_hypercube")
    
    SI = SourceIteration(data)
    SI.max_iter = 25
    SI.Run()

    

    