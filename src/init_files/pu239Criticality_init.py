#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon Apr  4 17:52:36 2022

@author: sampasmann
"""

import numpy as np
from src.functions.material import Material
from src.functions.mesh import Mesh
from src.init_files.reeds_solution import reeds_sol

class pu239Criticality_init:
    def __init__(self, N=2**10, Nx=100, generator="halton"):
        self.N = N
        self.Nx = Nx
        self.generator = generator
        self.totalDim = 2
        self.RB = 2*1.853722
        self.LB = 0#-1.853722
        self.G = 1
        self.right = False
        self.left = False
        self.material_code = "pu239_data"
        self.geometry = "slab"
        self.avg_scalar_flux = True
        self.edge_scalar_flux = False
        self.avg_angular_flux = False
        self.avg_current = False
        self.edge_current = False
        self.shannon_entropy = False
        self.save_data = True
        self.true_flux = np.array((False))
        self.mesh = Mesh(self.LB, self.RB, self.Nx)
        self.material = Material(self.material_code, self.geometry, self.mesh)
        self.source = np.ones((self.Nx,self.G))