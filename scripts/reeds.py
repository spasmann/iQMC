# -*- coding: utf-8 -*-
import sys
sys.path.append("../")
from src.init_files.reeds_init import ReedsInit
from src.solvers.solvers import LGMRES
import matplotlib.pyplot as plt

if __name__ == "__main__":
    N = 2**11
    Nx = 16*5
    generator = "halton"
    
    data = ReedsInit(N=N, Nx=Nx, generator=generator)
    phi = LGMRES(data)
    sol = data.true_flux
    
    plt.figure()
    plt.plot(range(Nx), phi)
    plt.plot(range(Nx), sol)

    
        
        