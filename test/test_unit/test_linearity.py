# -*- coding: utf-8 -*-
"""
Created on Wed Oct 19 10:43:20 2022

@author: Sam Pasmann
"""
import sys
sys.path.append("../../")
from src.input_files.reeds_init import ReedsInit
from src.solvers.fixed_source.maps import RHS
from src.functions.tallies import Tallies
from src.functions.sweep import Sweep
from src.functions.source import GetSource
import numpy as np


def fixed_source_sweep(phi_in, data):
    tallies         = Tallies(data)
    tallies.phi_avg = phi_in
    sweep           = Sweep(data)
    source          = GetSource(phi_in, data)    
    sweep.Run(tallies, source)
    phi_out         = tallies.phi_avg
    return phi_out

def fixed_source_mat_mul(phi_in, data, b):
    mmul = fixed_source_sweep(phi_in, data) - b
    return mmul

def test_fixed_source_linearity():
    data    = ReedsInit(N=2**7, Nx=18, generator="halton")
    G       = data.G
    Nx      = data.Nx
    b       = RHS(data)
    phi0    = np.ones((Nx,G))
    phi1    = fixed_source_sweep(phi0, data)
    r0      = phi1 - phi0
    v1      = phi0 + r0
    mv1     = fixed_source_mat_mul(v1, data, b)
    mphi0   = fixed_source_mat_mul(phi0, data, b)
    mr0     = fixed_source_mat_mul(r0, data, b)
    linerr  = mv1 - (mphi0 + mr0)
    assert(np.isclose(np.linalg.norm(linerr), 0.0))
    
if (__name__ == "__main__"):
    test_fixed_source_linearity()
    