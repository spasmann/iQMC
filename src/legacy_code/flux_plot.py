#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon Jun 28 14:09:33 2021

@author: sampasmann
"""
import numpy as np
import matplotlib.pyplot as plt
from cycles import cycles
from kmatrix import kmatrix

Geo = 1
#number of particles per cell
N = 5000
#number of cells
Nr = 10


#Nuclear data
nu = np.array([(2.84)])
nu = np.reshape(nu,(1,1))

Sig_f = np.array([(0.0816)])
Sig_f = np.reshape(Sig_f,(1,1))

Sig_c = np.array([(0.019584)])
Sig_c = np.reshape(Sig_c,(1,1))

Sig_s = np.array([(0.225216)])
Sig_s = np.reshape(Sig_s,(1,1,1))

Sig_t = np.array([(0.32640)])
Sig_t = np.reshape(Sig_t,(1,1))

X = np.array([(1.0)])
X = np.reshape(X,(1,1))

#critical radius of geo
if (Geo == 1):
    R = np.array([2.256751]) #slab
elif (Geo == 2):
    R = np.array([4.27996]) #cylinder
else:
    R = np.array([6.082547]) #sphere
    
#k_matrix, scalar_flux_m,H = kmatrix(Geo,N,R,Nr,Sig_t,Sig_f,Sig_s,Sig_c,X,nu)


k, scalar_flux,entropy = cycles(Geo,N,Nr,R, Sig_t,Sig_f,Sig_s,Sig_c,X,nu,weights = [],
                       inactive_cycles = 5, active_cycles = 20)



#flux plotting
plt.figure(dpi=200)
dr = R[0]/Nr
radii = np.linspace(dr/2,R[0]-dr/2,num=Nr)
plt.plot(radii,scalar_flux[0])
flux_r = np.array([0.25,0.5,0.75,1.0])*R[0]
if (Geo == 3):
    flux = np.array([0.93538006, 0.75575352, 0.49884364, 0.19222603])
elif (Geo ==2):
    flux = np.array([0.8093,0.2926])
    flux_r = np.array([0.5,1.0])*R[0]
else:
    flux = np.array([0.9701734, 0.8810540, 0.7318131, 0.4902592])

#for z in range(len(flux_r)):
plt.plot(flux_r,flux,'ro',label='test data')