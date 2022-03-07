#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon Mar  7 12:00:10 2022

@author: sampasmann
"""
import h5py
import numpy as np

def SaveData(init_data, tallies, fname = "default", path = "../saved_data/"):
    if (fname == "default"):
        fname = "{}-{}-{}-{}".format(init_data.material_code,
                                     init_data.generator,
                                     init_data.N,
                                     init_data.Nx)
    
    with h5py.File(path+fname, 'w') as f:
        f.create_dataset('N', data = init_data.N)
        f.create_dataset('Nx', data = init_data.Nx)
        f.create_dataset('phi_avg', data = tallies.phi_avg)
    
    return