#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon Apr  4 15:04:59 2022

@author: sampasmann
"""

import numpy as np
import os

def SHEM_361_data(Nx):
    script_dir = os.path.dirname(__file__) #<-- absolute dir the script is in
    rel_path = "../materials/SHEM-361.npz"
    abs_file_path = os.path.join(script_dir, rel_path)
    
    with np.load(abs_file_path) as data:
        data    = np.load('SHEM-361.npz')
        sigt    = data['SigmaT']
        sigs    = data['SigmaS']
        siga    = data['SigmaA']
        sigf    = data['SigmaF']
        nu      = data['nuSigmaF']
        chi     = data['chi_p']
        
    # Add leakage to make it subcritical
    sigt += 0.0036 # k = 0.99112    
    G = len(SigmaT) # number of groups

    # repeat xs for shape (Nx, G)
    sigs = np.tile(sigs,(Nx,1,1)) # (Nx, G, G)
    siga = np.tile(siga, (Nx,1))
    sigt = np.tile(sigt, (Nx,1))  
    sigf = np.tile(sigf, (Nx,1))
    nu = np.tile(nu,(Nx,G))
    chi = np.tile(chi,(Nx,G))
    
    return sigt, sigs, siga, sigf, chi, nu, G

if (__name__ == "__main__"):
    shem_361_data(20)