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

class PUa_1_0_SP_init:
    def __init__(self, N=2**10, Nx=100, generator="halton", source_tilt=False):
        np.random.seed(123456)
        self.keff               = 1.0
        self.N                  = N
        self.Nx                 = Nx
        self.generator          = generator
        self.source_tilt        = source_tilt
        self.totalDim           = 3
        self.RB                 = 6.082547
        self.LB                 = 0.0
        self.rng_seed           = 12345
        self.material_code      = "PUa_1_0"
        self.geometry           = "sphere"
        self.mode               = "eigenvalue"
        self.flux               = True
        self.save_data          = False
        self.right              = False
        self.left               = False
        self.RQMC               = False
        self.true_flux          = np.array((False))
        self.mesh               = Mesh(self.LB, self.RB, self.Nx)
        self.material           = Material(self.material_code, self.geometry, self.mesh)
        self.G                  = self.material.G
        self.fixed_source       = np.zeros((self.Nx,self.G))
        self.tallies            = Tallies(self)
        self.Nt                 = int(self.Nx*self.G)
        if (self.source_tilt):
            self.Nt = int(self.Nt*2)