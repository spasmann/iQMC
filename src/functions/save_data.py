#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon Mar  7 12:00:10 2022

@author: sampasmann
"""
import h5py
import numpy as np

def SaveData(init_data, SI, fname = "default", path = "../saved_data/"):
    if (fname == "default"):
        fname = "{}-{}-{}-{}".format(init_data.material_code,
                                     init_data.generator,
                                     init_data.N,
                                     init_data.Nx)
    
    with h5py.File(path+fname, 'w') as f:
        f.create_dataset('N', data = init_data.N)
        f.create_dataset('Nx', data = init_data.Nx)
        f.create_dataset('LB', data = init_data.LB)
        f.create_dataset('RB', data = init_data.RB)
        f.create_dataset('generator', data = init_data.generator)
        f.create_dataset('itt', data = SI.itt)
        f.create_dataset('phi_avg', data = SI.tallies.phi_avg)
        f.create_dataset('delta_flux', data = SI.norm_hist)
        if (init_data.true_flux.any()):
            f.create_dataset('true_flux', data = init_data.true_flux)
            f.create_dataset('error', data = SI.error)
    print("---------------------------------------")
    print("Simulation Data Saved at: ",path+fname)
    print("---------------------------------------")
    return