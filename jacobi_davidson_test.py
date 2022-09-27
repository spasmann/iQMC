#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Sep 22 15:33:47 2022

@author: sampasmann
"""

import numpy as np
import scipy.linalg as sp
import time
import sys
from scipy.sparse.linalg import lgmres, LinearOperator
from src.solvers.eigenvalue.maps import MatVec_data, MatVec
from src.input_files.PUa_1_0_SL_init import PUa_1_0_SL_init

Nx = 50
N = 1000
G = 1
Nv = int(Nx*G)
generator = "halton"
qmc_data = PUa_1_0_SL_init(N=N, Nx=Nx, generator=generator)
matvec_data = MatVec_data(qmc_data)
A           = LinearOperator((Nv,Nv), 
                              matvec=MatVec,
                              rmatvec=MatVec,
                              matmat= MatVec,
                              rmatmat=MatVec,
                              dtype=float)
b           = matvec_data[0]
phi0        = qmc_data.source
phi0        = np.reshape(phi0,(Nv,1))


# set up subspace "trial vectors"
k   = 8                  # number of initial guess vectors  
eig = 1                  # number of eignvalues to solve  
t   = np.eye(n,n)        # set of k unit vectors as guess
t   = t[:,(n-k):]  
V   = np.zeros((n,n))    # array of zeros to hold guess vec  
I   = np.eye(n)          # identity matrix same dimensions as A 
count = 1