#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon Apr  4 18:29:16 2022

@author: sampasmann
"""

import sys
sys.path.append("../")
from src.init_files.u235H2O_init import u235H2O_init
from src.functions.power_iteration import PowerIteration

if __name__ == "__main__":
    # initialize problem data
    N = 2**10
    data = u235H2O_init(N=2**10, Nx=10, generator="halton")

    SI = PowerIteration(data)
    SI.max_iter = 20
    # run source iteration
    SI.Run()