#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon Apr  4 18:29:16 2022

@author: sampasmann
"""

import sys
sys.path.append("../")
from src.init_files.URRb_H2Oa5_2_0_SL_init import URRb_H2Oa5_2_0_SL_init
from src.functions.power_iteration import PowerIteration

if __name__ == "__main__":
    # initialize problem data
    data = URRb_H2Oa5_2_0_SL_init(N=1000, Nx=20, generator="halton")
    data.save_data = False

    SI = PowerIteration(data)
    SI.max_iter = 20
    # run source iteration
    SI.Run()