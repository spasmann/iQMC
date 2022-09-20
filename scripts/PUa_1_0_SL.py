#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon Apr  4 17:52:16 2022

@author: sampasmann
"""

import sys
sys.path.append("../")
from src.input_files.PUa_1_0_SL_init import PUa_1_0_SL_init
from src.solvers.eigenvalue.solvers import PowerIteration

if __name__ == "__main__":
    
    # initialize problem data
    Nx = 20
    data = PUa_1_0_SL_init(N=10000, Nx=Nx, generator="random")
    data.save_data = False
    data.RQMC = False
    PI = PowerIteration(data)
    
    
