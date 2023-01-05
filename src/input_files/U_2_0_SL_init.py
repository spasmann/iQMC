#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon Apr  4 18:23:32 2022

@author: sampasmann
"""

import numpy as np
from src.functions.material import Material
from src.functions.mesh import Mesh
from src.functions.tallies import Tallies

class U_2_0_SL_init:
    def __init__(self, N=2**10, Nx=20, generator="halton",source_tilt=False):
        self.keff               = 1.0
        self.N                  = N
        self.Nx                 = Nx
        self.generator          = generator
        self.source_tilt        = source_tilt
        self.totalDim           = 2
        self.RB                 = 3.006375
        self.LB                 = -3.006375
        self.rng_seed           = 12345
        self.right              = False
        self.left               = False
        self.material_code      = "U_2_0"
        self.geometry           = "slab"
        self.mode               = "eigenvalue"
        self.flux               = True
        self.save_data          = False
        self.right              = False
        self.left               = False
        self.RQMC               = False
        self.true_flux          = np.array((False))
        self.mesh               = Mesh(self.LB, self.RB, self.Nx, self.geometry)
        self.material           = Material(self.material_code, self.geometry, self.mesh)
        self.G                  = self.material.G
        self.fixed_source       = np.zeros((self.Nx,self.G))
        self.tallies            = Tallies(self)
        self.Nt                 = int(self.Nx*self.G)
        if (self.source_tilt):
            self.Nt = int(self.Nt*2)
        