# -*- coding: utf-8 -*-
import sys, os
sys.path.append(os.getcwd()+"/../../../")
from src.input_files.reeds_init import ReedsInit
from src.solvers.fixed_source.solvers import FixedSource
import numpy as np

def test_reeds():
    """
    Test we are getting the same answers for using both a Krylov and Picard 
    method.
    
    Tests multimedia capabilites for Fixed-Source problems.
    """
    a           = np.genfromtxt("answers_gmres.csv",  delimiter=',')
    b           = np.genfromtxt("answers_picard.csv",  delimiter=',')
    N           = 2**6
    Nx          = 16
    maxit       = 10
    a           = np.reshape(a, (Nx,1))
    b           = np.reshape(b, (Nx,1))
    generator   = "sobol"
    solver      = "GMRES"
    data        = ReedsInit(N=N, Nx=Nx, generator=generator)
    c           = FixedSource(data,solver=solver,maxit=maxit,
                              report_progress = False)
    
    solver      = "Picard"
    data        = ReedsInit(N=N, Nx=Nx, generator=generator)
    d           = FixedSource(data,solver=solver,maxit=maxit,
                              report_progress = False)
    
    assert(np.allclose(a,c))
    assert(np.allclose(b,d))
    
    
if (__name__ == "__main__"):
    test_reeds()