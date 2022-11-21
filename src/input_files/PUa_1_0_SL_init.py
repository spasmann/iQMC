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

class PUa_1_0_SL_init:
    def __init__(self, N=2**10, Nx=100, generator="halton"):
        self.keff               = 1.0
        self.N                  = N
        self.Nx                 = Nx
        self.generator          = generator
        self.totalDim           = 2
        self.RB                 = 2.256751
        self.LB                 = -2.256751
        self.rng_seed           = 123456
        self.material_code      = "PUa_1_0"
        self.geometry           = "slab"
        self.mode               = "eigenvalue"
        self.flux               = True
        self.flux_derivative    = False
        self.source_tilt        = False
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