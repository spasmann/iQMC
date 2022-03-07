# -*- coding: utf-8 -*-

import sys
sys.path.append("../")
from src.init_files.garcia_init import GarciaInit
from src.functions.source_iteration import SourceIteration

import matplotlib.pyplot as plt


# initialize problem data
N = 2**10
data1 = GarciaInit(N=N, generator="random")
data1.material.c = 1.0
data2 = GarciaInit(N=N, generator="sobol")
data2.material.c = 1.0

# initialize source iteration
SI = SourceIteration(data1)
SI.max_iter = 20
# run source iteration
SI.Run()
phi1 = SI.tallies.phi_avg

SI = SourceIteration(data2)
SI.max_iter = 20
# run source iteration
SI.Run()
phi2 = SI.tallies.phi_avg
Nx = range(data2.Nx)

plt.plot(Nx, phi1, label='random')
plt.plot(Nx, phi2, label='sobol')

