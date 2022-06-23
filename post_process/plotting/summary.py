#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon May  9 14:57:15 2022

@author: sampasmann
"""
import matplotlib.pyplot as plt
import matplotlib
import numpy as np
import h5py

def summary(problem, generator, N, Nx):
    
    N = str(N)
    Nx = str(Nx)
    
    try:
        path = "../saved_data/"
        file = problem+"-"+generator+"-"+N+"-"+Nx
        f = h5py.File(path+file, 'r')    
    except:
        print("File DNE: ",path+file)
        return

    phi = f['phi_avg'][...]
    true = f['true_flux'][...]
    error = f['error'][...]
    
    MaxErr = abs((phi-true)).max()
    MaxRelErr = abs((phi-true)/true).max()
    MinErr = abs((phi-true)).min()
    MinRelErr = abs((phi-true)/true).min()
    
    print()
    print("----------------------------------------------------")
    print("          SUMMARY of "+file+"          ")
    print("----------------------------------------------------")
    print()
    print("Variables in H5PY File: ", list(f.keys()))
    print("- - - - - - - - - - - - - - - - - - - - - - - - - -")
    print("Maximum Error: ", MaxErr)
    print("Minimum Error: ", MinErr)
    print("- - - - - - - - - - - - - - - - - - - - - - - - - -")
    print("Maximum Relative Error: ",MaxRelErr)
    print("Minimum Relative Error: ",MinRelErr)
    print("- - - - - - - - - - - - - - - - - - - - - - - - - -")
    print("Max Error Stored in File: ",error.max())
    print("Min Error Stored in File: ",error.min())
    print("- - - - - - - - - - - - - - - - - - - - - - - - - -")
    print("Number of iterations: ", f['itt'][...])
    print("----------------------------------------------------")
    print("                   END OF SUMMARY                   ")
    print("----------------------------------------------------")

    
    return
