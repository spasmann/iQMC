#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon Apr  4 15:04:59 2022

@author: sampasmann
"""

import numpy as np
import os

def hdpe_data(material_code, Nx):
    script_dir = os.path.dirname(__file__) #<-- absolute dir the script is in
    rel_path = "../materials/HDPE/"
    abs_file_path = os.path.join(script_dir, rel_path)
    G = material_code
    D = np.genfromtxt(abs_file_path+"D_{}G_HDPE.csv".format(G), delimiter=",")
    siga = np.genfromtxt(abs_file_path+"Siga_{}G_HDPE.csv".format(G), delimiter=",")
    sigs = np.genfromtxt(abs_file_path+"Scat_{}G_HDPE.csv".format(G), delimiter=",")
    sigs = np.flip(sigs,1)
    sigs = np.flip(sigs)
    sigs = np.tile(sigs,(Nx,1,1))
    siga = np.tile(siga, (Nx,1))
    sigt = 1/(3*D)
    sigt = np.flip(np.tile(sigt, (Nx,1)))  # repeat sigt so its shape (Nx, G)
    
    return sigt, sigs, siga, G
