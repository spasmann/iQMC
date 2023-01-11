# -*- coding: utf-8 -*-
import sys
import os 
script_dir = os.path.dirname(__file__) #<-- absolute dir the script is in
rel_path = "../../../"
abs_file_path = os.path.join(script_dir, rel_path)
sys.path.append(abs_file_path)
from src.input_files.mg_init import MultiGroupInit
from src.solvers.fixed_source.solvers import FixedSource
import numpy as np

def test_hdpe_inf_slab():
    """
    Test we are getting the same answers for using both a Krylov and Picard 
    method.
    
    Tests multigroup capabilites for Fixed-Source problems.
    """
    a           = np.genfromtxt("./test/test_integration/test_hdpe_inf_slab/answers_gmres.csv",  delimiter=',')
    b           = np.genfromtxt("./test/test_integration/test_hdpe_inf_slab/answers_picard.csv",  delimiter=',')
    N           = 2**6
    Nx          = 10
    G           = 12
    maxit       = 10
    a           = np.reshape(a, (Nx,G))
    b           = np.reshape(b, (Nx,G))
    generator   = "sobol"
    solver      = "GMRES"
    data        = MultiGroupInit(numGroups=G, N=N, Nx=Nx, generator=generator)
    c           = FixedSource(data,solver=solver,maxit=maxit,
                              report_progress = False)
    
    solver      = "Picard"
    data        = MultiGroupInit(numGroups=G, N=N, Nx=Nx, generator=generator)
    d           = FixedSource(data,solver=solver,maxit=maxit,
                              report_progress = False)
    
    assert(np.allclose(a,c))
    assert(np.allclose(b,d))
    
    
if (__name__ == "__main__"):
    test_hdpe_inf_slab()