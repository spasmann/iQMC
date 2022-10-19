# -*- coding: utf-8 -*-
"""
Created on Wed Oct 19 10:43:20 2022

@author: Sam Pasmann
"""
import sys
sys.path.append("../../")
from src.input_files.reeds_init import ReedsInit
from src.input_files.Ua_1_0_SL_init import Ua_1_0_SL_init
from src.solvers.fixed_source.maps import RHS as fixed_source_RHS
from src.solvers.eigenvalue.davidson import PreConditioner
from src.functions.tallies import Tallies
from src.functions.sweep import Sweep
from src.functions.source import GetSource, GetCriticalitySource
import numpy as np
import scipy.linalg as sp

from src.solvers.eigenvalue.davidson import AxV, BxV

def fixed_source_sweep(phi_in, data):
    tallies         = Tallies(data)
    #tallies.phi_avg = phi_in
    sweep           = Sweep(data)
    source          = GetSource(phi_in, data)    
    sweep.Run(tallies, source)
    phi_out         = tallies.phi_avg
    return phi_out

def fixed_source_mat_mul(phi_in, data, b):
    mmul = fixed_source_sweep(phi_in, data) - b
    return mmul

def eigenval_sweep(phi_f_in, phi_s_in, data):
    tallies     = Tallies(data)
    # tallies.phi_avg ?
    data.phi_f  = phi_f_in
    source      = GetCriticalitySource(phi_f_in, phi_s_in, data)
    sweep       = Sweep(data) 
    sweep.Run(tallies, source) # QMC sweep
    phi_out     = tallies.phi_avg
    return phi_out

def eigenval_mat_mul(phi_in, data):
    axv      = AxV(phi_in, data)
    bxv      = BxV(phi_in, data)
    AV       = np.dot(phi_in.T, axv) # Scattering linear operator
    BV       = np.dot(phi_in.T, bxv) # Fission linear operator
    [Lambda, w] = sp.eig(AV, b=BV)
    Lambda   = Lambda.real
    mmul = AV - Lambda*BV
    return mmul

def test_fixed_source_linearity():
    """
    Generally assuring that (Ax - b) = 0
    """
    data    = ReedsInit(N=2**7, Nx=18, generator="halton")
    G       = data.G
    Nx      = data.Nx
    b       = fixed_source_RHS(data)
    phi0    = np.ones((Nx,G))
    phi1    = fixed_source_sweep(phi0, data)
    r0      = phi1 - phi0
    v1      = phi0 + r0
    mv1     = fixed_source_mat_mul(v1, data, b)
    mphi0   = fixed_source_mat_mul(phi0, data, b)
    mr0     = fixed_source_mat_mul(r0, data, b)
    linerr  = mv1 - (mphi0 + mr0)
    assert(np.isclose(np.linalg.norm(linerr), 0.0))

def test_eigenval_linearity():
    """
    Generally assuring that (Ax - Lambda*Bx) = 0
    
    This tests functions for both PI and Davidson.
    """
    data        = Ua_1_0_SL_init(N=2**7, Nx=10, generator="sobol")
    G           = data.G
    Nx          = data.Nx
    phi0        = np.ones((Nx,G))
    data.phi_f  = phi0
    phi1        = eigenval_sweep(phi0, phi0, data)
    r0          = phi1 - phi0
    v1          = phi0 + r0
    mv1         = eigenval_mat_mul(v1, data)
    mphi0       = eigenval_mat_mul(phi0, data)
    mr0         = eigenval_mat_mul(r0, data)
    linerr      = mv1 - (mphi0 + mr0)
    assert(np.isclose(np.linalg.norm(linerr), 0.0))
    
#def test_davidson_preconditioner_liearity():
    # should match structure of fixed_source_linearity
    # since it is a series of scattering source sweeps

if (__name__ == "__main__"):
    test_fixed_source_linearity()
    test_eigenval_linearity()
    