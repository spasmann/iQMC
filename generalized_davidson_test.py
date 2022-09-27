#!/usr/bin/env python3
# -*- coding: utf-8 -*-

import numpy as np
import scipy.linalg as sp
import time
import sys
sys.path.append("./iQMC/")
from src.solvers.eigenvalue.maps import SI_Map
from src.input_files.PUa_1_0_SL_init import PUa_1_0_SL_init
from src.solvers.fixed_source.solvers import FixedSource

import matplotlib.pyplot as plt

def AxV(V, qmc_data):
    """
    Linear operator for scattering term (I-L^(-1)S)*phi
    """
    Nx  = qmc_data.Nx
    G   = qmc_data.G
    Nv  = Nx*G
    zed = np.zeros((Nx,G))
    axv = np.empty(V.shape)
    for i in range(V.shape[1]):
        phi_in = np.reshape(V[:,i], (Nv,1))
        axv[:,i] = phi_in[:,0] - SI_Map(zed, phi_in, qmc_data)[:,0]
        
    return axv

def BxV(V, qmc_data):
    """
    Linear operator for fission term (L^(-1)F*phi)
    """
    Nx  = qmc_data.Nx
    G   = qmc_data.G
    Nv  = Nx*G
    zed = np.zeros((Nx,G))
    bxv = np.empty(V.shape)
    for i in range(V.shape[1]):
        phi_in = np.reshape(V[:,i], (Nv,1))
        bxv[:,i]  = SI_Map(phi_in, zed, qmc_data)[:,0]
        
    return bxv

def PreConditioner(V,qmc_data, solver="GMRES"):
    """
    Linear operator approximation of L^(-1)S
    """
    maxit   = 10
    tol     = 1e-6
    Nx      = qmc_data.Nx
    G       = qmc_data.G
    Nv      = Nx*G
    t       = np.empty(V.shape)
    for i in range(V.shape[1]):
        qmc_data.source = np.reshape(V[:,i], (Nv,1))
        t[:,i] = FixedSource(qmc_data,solver=solver, maxit=maxit, tol=tol, report_progress=False, save_data=False)[:,0]
    return t

def Residual(V, Lambda, qmc_data):
    
    return r

def Gram(V,u):
    w1  = u - np.dot(V,np.dot(V.T,u))
    v1  = w1 / np.linalg.norm(w1)
    w2  = v1 - np.dot(V,np.dot(V.T,v1))
    v2  = w2 / np.linalg.norm(w2)
    V   = np.append(V, v2, axis=1)
    return V

# Problem setup
Nx          = 20
N           = 2**10
G           = 1
Nv          = int(Nx*G)
generator   = "halton"
qmc_data    = PUa_1_0_SL_init(N=N, Nx=Nx, generator=generator)
phi0        = qmc_data.source
phi0        = np.reshape(phi0,(Nv,1))

# Davidson Parameters
start = time.time()
u       = phi0.copy()
V       = np.array(u/np.linalg.norm(u).T) # orthonormalize initial guess
Lambda0 = 1.0
k_old   = 0.0
dk      = 1.0
r0      = AxV(V, qmc_data) - Lambda0*BxV(V, qmc_data) # (A - keff*B)V_0
r       = r0
tol     = 1e-9
itt     = 1
maxit   = 30
l       = 1 # compute "l" largest eigenpairs
m       = 4 # restart parameter (ie maximum size of V)

# Davidson Routine
while (itt <= maxit) and (dk>=tol):
    print(" Davidson Iteration: ", itt)
    AV          = np.dot(V.T, AxV(V, qmc_data)) # Scattering linear operator
    BV          = np.dot(V.T, BxV(V, qmc_data)) # Fission linear operator
    [Lambda, w] = sp.eig(AV, b=BV)  # solve for eigenvalues and vectors
    idx         = Lambda.argsort()  # get indices of eigenvalues from smallest to largest
    Lambda      = Lambda[idx]       # sort eigenvalues from smalles to largest
    Lambda      = Lambda[:l].real   # take the real component of the l largest eigenvalues
    k           = 1/Lambda
    dk          = abs(k - k_old)
    print("dk: ",dk)
    print(V.shape)
    k_old       = k
    w           = w[:,idx]          # sort corresponding eigenvector
    w           = w[:,:l].real      # take the l largest eigenvectors
    u           = np.dot(V,w)       # Ritz vectors
    r           = AxV(u, qmc_data) - Lambda*BxV(u, qmc_data) # residual
    t           = PreConditioner(r, qmc_data)
    if (V.shape[1] <= m-l ):
        V = Gram(V,t) # appends new orthogonalization to V
    else:
        V = Gram(u,t) # "restarts" by appending to a new array 
    if (itt==maxit):
        print("Maximum number of iterations")
        break
    itt += 1

stop = time.time()
print("Davidson took: ", stop-start)
keff = 1/Lambda
phi  = V[:,0]

plt.plot(range(Nx),phi)
