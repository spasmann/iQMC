#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon Apr  4 18:23:32 2022

@author: sampasmann
"""

import numpy as np
from src.functions.material import Material
from src.functions.mesh import Mesh

class URRa_2_0_SL_init:
    def __init__(self, N=2**10, Nx=100, generator="halton"):
        np.random.seed(123456)
        self.keff               = 1.0
        self.N                  = N
        self.Nx                 = Nx
        self.generator          = generator
        self.totalDim           = 2
        self.RB                 = 7.566853 
        self.LB                 = 7.566853 
        self.right              = False
        self.left               = False
        self.material_code      = "URRa_2_0_SL"
        self.geometry           = "slab"
        self.avg_scalar_flux    = True
        self.edge_scalar_flux   = False
        self.avg_angular_flux   = False
        self.avg_current        = False
        self.edge_current       = False
        self.shannon_entropy    = False
        self.save_data          = True
        self.moment_match       = False
        self.true_flux          = np.array((False))
        self.mesh               = Mesh(self.LB, self.RB, self.Nx)
        self.material           = Material(self.material_code, self.geometry, self.mesh)
        self.G                  = self.material.G
        self.source             = np.ones((self.Nx,self.G))#np.random.random(size=(self.Nx,self.G))
        self.phi_f              = np.ones((self.Nx,self.G))#np.random.random(size=(self.Nx,self.G))