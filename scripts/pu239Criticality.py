#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon Apr  4 17:52:16 2022

@author: sampasmann
"""

import sys
sys.path.append("../")
from src.init_files.pu239Criticality_init import pu239Criticality_init
from src.functions.power_iteration import PowerIteration

if __name__ == "__main__":
    # initialize problem data
    data = pu239Criticality_init(N=10000, Nx=100, generator="halton")

    SI = PowerIteration(data)
    SI.max_iter = 20
    # run source iteration
    SI.Run()