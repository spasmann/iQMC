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
from src.solvers.eigenvalue.davidson import AxV, BxV
from src.solvers.eigenvalue.davidson import PreConditioner
from src.solvers.eigenvalue.maps import SI_Map
from src.solvers.eigenvalue.maps import RHS as eigenval_RHS
from src.solvers.eigenvalue.solvers import UpdateK
from src.functions.tallies import Tallies
from src.functions.sweep import Sweep
from src.functions.source import GetSource, GetCriticalitySource
import numpy as np
import scipy.linalg as sp


# =============================================================================
# Mapping Functions
# =============================================================================

def fixed_source_sweep(phi_in, data):
    tallies         = Tallies(data)
    sweep           = Sweep(data)
    source          = GetSource(phi_in, data)    
    sweep.Run(tallies, source)
    phi_out         = tallies.phi_avg
    return phi_out

def fixed_source_mat_mul(phi_in, data, b):
    res = phi_in - (fixed_source_sweep(phi_in, data) - b)
    return res

def eigenval_sweep(phi_f_in, phi_s_in, data):
    tallies     = Tallies(data)
    source      = GetCriticalitySource(phi_f_in, phi_s_in, data)
    sweep       = Sweep(data) 
    sweep.Run(tallies, source) # QMC sweep
    phi_out     = tallies.phi_avg
    return phi_out

def eigenval_mat_mul(phi_in, data):
    phi_in   = PreConditioner(phi_in, data, numSweeps=30)
    axv      = AxV(phi_in, data)
    bxv      = BxV(phi_in, data)
    AV       = np.dot(phi_in.T, axv) # Scattering linear operator
    BV       = np.dot(phi_in.T, bxv) # Fission linear operator
    [Lambda, w] = sp.eig(AV, b=BV)
    Lambda   = Lambda.real
    res = axv - Lambda*bxv
    return res

def eigenval_mat_mul_2(phi_in, data, b):
    res = phi_in - (eigenval_sweep(phi_in, phi_in, data) - b)
    return res

# =============================================================================
# Linearity Tests
# =============================================================================

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
    linerr  = np.linalg.norm(mv1 - (mphi0 + mr0))
    print(linerr)
    assert(np.isclose(linerr, 0.0))

def test_eigenval_linearity():
    """
    Generally assuring that (Ax - Lambda*Bx) = 0
    
    This tests functions for both PI and Davidson.
    """
    data        = Ua_1_0_SL_init(N=2**7, Nx=10, generator="sobol")
    G           = data.G
    Nx          = data.Nx
    phi0        = np.ones((Nx,G))
    phi1        = eigenval_sweep(phi0, phi0, data)
    r0          = phi1 - phi0
    v1          = phi0 + r0
    mv1         = eigenval_mat_mul(v1, data)
    mphi0       = eigenval_mat_mul(phi0, data)
    mr0         = eigenval_mat_mul(r0, data)
    linerr  = np.linalg.norm(mv1 - (mphi0 + mr0))
    print(linerr)
    assert(np.isclose(linerr, 0.0))
    
def test_davidson_preconditioner_liearity():
    data    = Ua_1_0_SL_init(N=2**7, Nx=10, generator="halton")
    sweeps  = 1
    G       = data.G
    Nx      = data.Nx
    phi0    = np.ones((Nx,G))
    phi1    = PreConditioner(phi0, data, numSweeps=sweeps)
    r0      = phi1 - phi0
    v1      = phi0 + r0
    mv1     = PreConditioner(v1, data, numSweeps=sweeps)
    mphi0   = PreConditioner(phi0, data, numSweeps=sweeps)
    mr0     = PreConditioner(r0, data, numSweeps=sweeps)
    linerr  = np.linalg.norm(mv1 - (mphi0 + mr0))
    print(linerr)
    assert(np.isclose(linerr, 0.0))

if (__name__ == "__main__"):
    test_fixed_source_linearity()
    test_eigenval_linearity()
    test_davidson_preconditioner_liearity()
    