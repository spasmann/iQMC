#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon Mar  7 12:00:10 2022

@author: sampasmann
"""
import h5py
import numpy as np
import os


def SaveData(qmc_data, sim_data, fname="default",
             path="../post_process/saved_data/"):
    if (fname == "default"):
        fname = "{}-{}-{}-{}-{}".format(qmc_data.material_code,
                                        qmc_data.generator,
                                        qmc_data.N,
                                        qmc_data.Nx,
                                        sim_data["nproc"])

    with h5py.File(path + fname, 'w') as f:
        f.create_dataset('N', data=qmc_data.N)
        f.create_dataset('Nx', data=qmc_data.Nx)
        f.create_dataset('LB', data=qmc_data.LB)
        f.create_dataset('RB', data=qmc_data.RB)
        f.create_dataset('generator', data=qmc_data.generator)
        f.create_dataset('phi_avg', data=sim_data["phi"])
        f.create_dataset('run_time', data=sim_data["run_time"])
        f.create_dataset('tolerance', data=sim_data["tolerance"])
        f.create_dataset('nproc', data=sim_data["nproc"])
        # f.create_dataset('delta_flux', data = SI.norm_hist)
        if (qmc_data.true_flux.any()):
            f.create_dataset('true_flux', data=qmc_data.true_flux)
    print("---------------------------------------")
    print("Simulation Data Saved at: ", os.getcwd() + "/" + path + fname)
    print("---------------------------------------")
    return
