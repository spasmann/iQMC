# -*- coding: utf-8 -*-

"""
Created on Mon May  9 14:57:15 2022

@author: sampasmann

In this script I try to dynamically increase the spatial cells when the
convergance rate is not within a speficied relative tolerance of the expected
1/N given a fixed number of particles. 

ie it will keep running the same number
of particles, increasing the spatial cells, until within the tolerance. If its 
within the tolerance on the first try, then increase the number of particles and 
go again.
"""

import sys
sys.path.append("../")
from src.init_files.reeds_init import ReedsInit
from src.functions.source_iteration import SourceIteration

import h5py
import math

if __name__ == "__main__":
    Nvals = [2**10,2**11,2**12,2**13,2**14,2**15]
    Nx = 128
    NList = []
    NxList = []
    count = 0
    N = Nvals[count]
    problem = "reeds_data"
    generator = "halton"
    path = "../saved_data/"
    rel_tol = 0.3 
    
    # fetch error from first problem
    name = problem+"-"+generator+"-"+str(N)+"-"+str(Nx)
    f = h5py.File(path+name,'r')
    error1 = f['error'][:][-1]
    f.close()
    count += 1
    # increase N
    N = Nvals[count]
    
    while (N <= Nvals[-1]):
        # run simulation
        data = ReedsInit(N=N, Nx=Nx, generator=generator)
        SI = SourceIteration(data)
        SI.max_iter = 100
        SI.Run()
        # read final error from simulation
        name = problem+"-"+generator+"-"+str(N)+"-"+str(Nx)
        f = h5py.File(path+name,'r')
        error2 = f['error'][:][-1]
        f.close()
        errortheo = error1*N/Nvals[count+1] # theoretical error
        print("---------------------------------")
        print(name, ": Error = ",error2)
        print("Theoretical Error = ", errortheo)
        print("---------------------------------")
        if (math.isclose(error2,errortheo,rel_tol=rel_tol)):
            print("Within Convergance Tolerance, increasing N:")
            print(N,"->",Nvals[count+1])
            error1 = error2
            NList.append(N)
            NxList.append(Nx)
            count += 1
            N = Nvals[count]
        else:
            print("Outside Convergance Tolerance, increasing Nx:")
            print(Nx,"->",Nx+16)
            Nx += 16

        if (error2 > error1):
            break

    