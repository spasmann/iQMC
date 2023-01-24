#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon Apr  4 17:52:36 2022

@author: sampasmann
"""

import numpy as np
from src.functions.material import Material
from src.functions.mesh import Mesh
from src.functions.tallies import Tallies

class Ua_1_0_SL_init:
    def __init__(self, N=2**10, Nx=100, generator="halton",source_tilt=False,
                 seed=12345, RQMC=False):
        self.keff               = 1.0
        self.N                  = N
        self.Nx                 = Nx
        self.generator          = generator
        self.source_tilt        = source_tilt
        self.rng_seed           = seed
        self.RQMC               = RQMC
        self.totalDim           = 2
        self.RB                 = 2.872934 
        self.LB                 = -2.872934 
        self.right              = False
        self.left               = False
        self.material_code      = "Ua_1_0"
        self.geometry           = "slab"
        self.mode               = "eigenvalue"
        self.flux               = True
        self.save_data          = False
        self.right              = False
        self.left               = False
        self.true_flux          = np.array((False))
        self.mesh               = Mesh(self.LB, self.RB, self.Nx, self.geometry)
        self.material           = Material(self.material_code, self.geometry, self.mesh)
        self.G                  = self.material.G
        self.fixed_source       = np.zeros((self.Nx,self.G))
        self.FixedSource        = lambda x,cell: self.fixed_source[cell,:]
        self.tallies            = Tallies(self)
        self.Nt                 = int(self.Nx*self.G)
        if (self.source_tilt):
            self.Nt = int(self.Nt*2)