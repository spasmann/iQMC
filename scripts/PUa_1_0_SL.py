#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon Apr  4 17:52:16 2022

@author: sampasmann
"""

import sys
sys.path.append("../")
from src.init_files.PUa_1_0_SL_init import PUa_1_0_SL_init
from src.functions.power_iteration import PowerIteration

if __name__ == "__main__":
    
    # initialize problem data
    data = PUa_1_0_SL_init(N=1000, Nx=50, generator="halton")
    data.save_data = False
    data.RQMC = False

    PI = PowerIteration(data)
    PI.max_iter = 10
    # run source iteration
    PI.Run()
    
    
