# -*- coding: utf-8 -*-
import sys
sys.path.append("../")
from src.init_files.reeds_init import ReedsInit
from src.functions.source_iteration import SourceIteration

if __name__ == "__main__":
    N = 2**11
    Nx = 40
    generator = "halton"
    data = ReedsInit(N=N, Nx=Nx, generator=generator, LB=0.0)
    
    SI = SourceIteration(data)
    SI.max_iter = 100
    SI.Run()
    
        
        