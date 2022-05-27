#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon Apr  4 18:29:16 2022

@author: sampasmann
"""

import sys
sys.path.append("../")
from src.init_files.PUa_H2O_1_0_SL_init import PUa_H2O_1_0_SL_init
from src.functions.power_iteration import PowerIteration

if __name__ == "__main__":
    # initialize problem data
    data = PUa_H2O_1_0_SL_init(N=1000, Nx=25, generator="halton")
    data.save_data = False

    SI = PowerIteration(data)
    SI.max_iter = 10
    # run source iteration
    SI.Run()