# -*- coding: utf-8 -*-
"""
Created on Wed Oct 19 10:43:20 2022

@author: Sam Pasmann
"""
import sys
import os 
script_dir = os.path.dirname(__file__) #<-- absolute dir the script is in
rel_path = "../../"
abs_file_path = os.path.join(script_dir, rel_path)
sys.path.append(abs_file_path)
import numpy as np
from src.functions.sweep import Sweep
from src.functions.tallies import Tallies
from src.input_files.reeds_init import ReedsInit
from src.input_files.Ua_1_0_SL_init import Ua_1_0_SL_init
from src.solvers.eigenvalue.davidson import AxV, BxV
from src.solvers.eigenvalue.davidson import PreConditioner
from src.solvers.fixed_source.maps import RHS as fixed_source_RHS
from src.functions.source import GetSource

# =============================================================================
# Mapping Functions
# =============================================================================

def fixed_source_sweep(phi_in, data):
    data.tallies.q  = GetSource(phi_in, data)
    sweep           = Sweep(data)   
    sweep.Run(data)
    phi_out         = data.tallies.phi_avg
    return phi_out

def fixed_source_mat_mul(phi_in, data, b):
    res = (fixed_source_sweep(phi_in, data) - b)
    return res

def eigenval_sweep(phi_f_in, phi_s_in, data):
    data.tallies.q  = GetSource(phi_s_in, data, phi_avg_f=phi_f_in)
    sweep           = Sweep(data) 
    sweep.Run(data) # QMC sweep
    phi_out         = data.tallies.phi_avg
    return phi_out

def standard_eigenval_mat_mul(phi_f, phi_s, data, b):
    res = (eigenval_sweep(phi_f, phi_s, data) - b)
    return res

# =============================================================================
# Linearity Tests
# =============================================================================

def test_fixed_source_linearity():
    """
    Test linearity of the fixed source qmc sweep functions
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
    #print(linerr)
    assert(np.isclose(linerr, 0.0))
    
def test_standard_eigenval_linearity():
    data    = Ua_1_0_SL_init(N=2**7, Nx=18, generator="halton")
    G       = data.G
    Nx      = data.Nx
    zed     = np.zeros((Nx,G))
    b       = eigenval_sweep(zed, zed, data)
    phi0    = np.ones((Nx,G))
    phi1    = eigenval_sweep(phi0, phi0, data)
    r0      = phi1 - phi0
    v1      = phi0 + r0
    mv1     = standard_eigenval_mat_mul(v1, v1, data, b)
    mphi0   = standard_eigenval_mat_mul(phi0, phi0, data, b)
    mr0     = standard_eigenval_mat_mul(r0, r0, data, b)
    linerr  = np.linalg.norm(mv1 - (mphi0 + mr0))
    #print(linerr)
    assert(np.isclose(linerr, 0.0))
        

def test_davidson_preconditioner_liearity():
    data    = Ua_1_0_SL_init(N=2**7, Nx=10, generator="halton")
    sweeps  = 8
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
    #print(linerr)
    assert(np.isclose(linerr, 0.0))


if (__name__ == "__main__"):
    test_fixed_source_linearity()
    test_standard_eigenval_linearity()
    test_davidson_preconditioner_liearity()
    