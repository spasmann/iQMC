# -*- coding: utf-8 -*-
"""
Created on Thu Jan 12 16:22:00 2023

@author: Sam
"""
import sys
sys.path.append("../")
from src.input_files.PUa_1_0_SL_init import PUa_1_0_SL_init
from src.input_files.PUa_1_0_SP_init import PUa_1_0_SP_init
from src.solvers.eigenvalue.solvers import PowerIteration, Davidson
import matplotlib.pyplot as plt
import numpy as np

if __name__ == "__main__":
    # N_list      = np.array((10,50,100,500,1000,5000))
    N_list      = np.array((250,500,2500,5000))
    # N_list      = 10*np.array((1,2,3,4,5,6))
    J           = 1
    Nx          = 12
    maxit       = 50
    tol         = 1e-3
    solver      = "LGMRES"
    seeds       = np.array((1000000*np.random.random(J)),dtype=int)
    k1          = np.zeros((N_list.size, J))
    k2          = np.zeros((N_list.size, J))
    
    avg1        = np.zeros(J)
    avg2        = np.zeros(J)
    
    for i in range(N_list.size):
        for j in range(J):
            N           = N_list[i]
            seed        = seeds[j]
            
            data1       = PUa_1_0_SP_init(N=N, Nx=Nx, generator="halton",
                                          source_tilt=True,
                                          RQMC = True,
                                          seed = seed)
            phi1, khist1, itt1 = PowerIteration(data1,
                                            solver=solver,
                                            max_outter_itt=maxit, 
                                            max_inner_itt=maxit, 
                                            outter_tol=tol,
                                            inner_tol=tol,
                                            report_progress=False)
            k1[i,j] = khist1[-1]

            
            data2       = PUa_1_0_SP_init(N=N, Nx=Nx, generator="random",
                                          source_tilt = True,
                                          RQMC = True,
                                          seed = seed)
            phi2, khist2, itt2 = PowerIteration(data2,
                                            solver=solver,
                                            max_outter_itt=maxit, 
                                            max_inner_itt=maxit, 
                                            outter_tol=tol,
                                            inner_tol=tol,
                                            report_progress=False)
            k2[i,j] = khist2[-1]

        

# =============================================================================
# basic plot
# =============================================================================

kerr1 = (abs(1.0-k1)/k1.sum()).sum(axis=1)
kerr2 = (abs(1.0-k2)/k2.sum()).sum(axis=1)
        
count1 = 0
count2 = -1

plt.figure(dpi=300,figsize=(8,5))
plt.title(r'$k_\mathrm{eff}$ Relative Error')

plt.plot(N_list[count1:], kerr2[count1]*np.sqrt(N_list[count1]/N_list[count1:]), 'b--',label='$N^{-1/2}$')
plt.plot(N_list[count1:], kerr2[count1:], 'b-o', label='MC')

plt.plot(N_list[count1:], kerr1[count1]*(N_list[count1]/N_list[count1:]), 'g--',label='$N^{-1}$')
plt.plot(N_list[count1:], kerr1[count1:], 'g-^', label='QMC')

plt.yscale('log')
plt.xscale('log')
plt.legend()

# =============================================================================
# plot source
# =============================================================================
# data = data2
# source_tilt = data.source_tilt

# plt.figure(figsize=(6,4),dpi=300)
# plt.title('Larsen Source')
# q    = data.tallies.q
# mesh = data.mesh.edges
# mid  = data.mesh.midpoints
# x    = np.linspace(data.LB, data.RB, num=1000)
# n    = len(x)
# conditions = [(mesh[i] <= x) & (x <= mesh[i+1]) for i in range(Nx)]
# y1 = np.piecewise(x, conditions, q)
# plt.plot(x,y1, label=r'$a_j$')
# if (source_tilt):
#     qdot        = data.tallies.qdot
#     y2          = np.zeros_like(x)
#     for i in range(n):
#         zone = data.mesh.GetZone([x[i],0,0], [0,0,0])
#         x_mid = data.mesh.midpoints[zone]
#         y2[i] = q[zone] + qdot[zone]*(x[i] - x_mid)
#     plt.plot(x,y2,label=r'$a_j + b_j(x)$')
# for i in range(len(mesh)):
#     plt.axvline(mesh[i],linestyle='-',color='black')
# plt.legend()
# plt.tight_layout()
        
    
# =============================================================================
# subplots grouped by MC/QMC
# =============================================================================

# plt.subplots(1,2, sharex=False, sharey=True, figsize=(9,4), dpi=300)
# plt.suptitle('Larsen et al. Infinite Linear Source Problem')
# count = 10

# plt.subplot(121)
# plt.title('QMC')
# plt.plot(N_list[:count], err1[0]*N_list[0]/N_list[:count], 'k-',label='$N^{-1}$')
# plt.plot(N_list[:count], err3[0]*N_list[0]/N_list[:count], 'k-')
# plt.plot(N_list[:count], err1[:count], 'b-o', label='Constant Source')
# plt.plot(N_list[:count], err3[:count], 'g--^', label='Linear Source')
# plt.yscale('log')
# plt.xscale('log')
# plt.ylabel(r'$\phi$ Relative Error')
# plt.xlabel('N')
# plt.legend()

# plt.subplot(122)
# plt.title('MC')
# plt.plot(N_list[:count], err2[0]*np.sqrt(N_list[0]/N_list[:count]), 'k-',label=r'$N^{-1/2}$')
# plt.plot(N_list[:count], err2[:count], 'b-o', label='Constant Source')
# plt.plot(N_list[:count], err4[:count], 'g--^', label='Linear Source')
# plt.yscale('log')
# plt.xscale('log')
# plt.xlabel('N')
# plt.legend()

# =============================================================================
# subplots grouped by Constant/Linear Source
# =============================================================================


# plt.subplots(1,2, sharex=False, sharey=True, figsize=(9,4), dpi=300)
# plt.suptitle('Larsen et al. Infinite Linear Source Problem')
# count = 10

# plt.subplot(121)
# plt.title('Constant Source')
# plt.plot(N_list[:count], err1[0]*N_list[0]/N_list[:count], 'k-',label='$N^{-1}$')
# plt.plot(N_list[:count], err2[0]*np.sqrt(N_list[0]/N_list[:count]), 'k--',label='$N^{-1/2}$')
# plt.plot(N_list[:count], err1[:count], 'b-o', label='QMC')
# plt.plot(N_list[:count], err2[:count], 'g--^', label='MC')
# plt.yscale('log')
# plt.xscale('log')
# plt.ylabel(r'$\phi$ Relative Error')
# plt.xlabel('N')
# plt.legend()

# plt.subplot(122)
# plt.title('Linear Source')
# plt.plot(N_list[:count], err3[0]*N_list[0]/N_list[:count], 'k-',label='$N^{-1}$')
# plt.plot(N_list[:count], err4[0]*np.sqrt(N_list[0]/N_list[:count]), 'k--',label='$N^{-1/2}$')
# plt.plot(N_list[:count], err3[:count], 'b-o', label='QMC')
# plt.plot(N_list[:count], err4[:count], 'g--^', label='MC')
# plt.yscale('log')
# plt.xscale('log')
# plt.xlabel('N')
# plt.legend()
