#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Apr  7 13:31:15 2022

@author: sampasmann
"""

import sys
sys.path.append("../")
from src.init_files.URRa_2_0_SL_init import URRa_2_0_SL_init
from src.functions.power_iteration import PowerIteration

if __name__ == "__main__":
    # initialize problem data
    data = URRa_2_0_SL_init(N=4000, Nx=int(100*7.566853), generator="halton")
    data.save_data = False

    SI = PowerIteration(data)
    SI.max_iter = 10
    # run source iteration
    SI.Run()