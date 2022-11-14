#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon Apr  4 17:52:36 2022

@author: sampasmann
"""

import numpy as np
from src.functions.material import Material
from src.functions.mesh import Mesh

class PUa_1_0_CY_init:
    def __init__(self, N=2**10, Nx=100, generator="halton"):
        np.random.seed(123456)
        self.keff               = 1.0
        self.N                  = N
        self.Nx                 = Nx
        self.generator          = generator
        self.totalDim           = 3
        self.RB                 = 4.279960
        self.LB                 = 0.0
        self.right              = False
        self.left               = False
        self.material_code      = "PUa_1_0"
        self.geometry           = "cylinder"
        self.flux               = True
        self.flux_derivative    = False
        self.save_data          = False
        self.moment_match       = False
        self.true_flux          = np.array((False))
        self.mesh               = Mesh(self.LB, self.RB, self.Nx)
        self.material           = Material(self.material_code, self.geometry, self.mesh)
        self.G                  = self.material.G
        self.source             = np.zeros((self.Nx,self.G))#np.random.random(size=(self.Nx,self.G))
        self.phi_f              = np.ones((self.Nx,self.G))#np.random.random(size=(self.Nx,self.G))