# -*- coding: utf-8 -*-

import numpy as np
from src.functions.material import Material
from src.functions.mesh import Mesh

class TestInit:
    def __init__(self, N=1000, Nx=20, generator="random"):
        self.N = N
        self.Nx = Nx
        self.generator = generator
        self.totalDim = 5
        self.RB = 5
        self.LB = 0
        self.G = 1
        self.right = False
        self.left = False
        self.material_code = "test_data"
        self.geometry = "slab"
        self.mesh = Mesh(self.Nx, np.array((self.RB,)))
        self.material = Material(self.material_code, self.geometry, self.mesh)
        self.source = np.ones((self.Nx,self.G))
        self.avg_scalar_flux = True
        self.edge_scalar_flux = False
        self.avg_angular_flux = False
        self.avg_current = False
        self.edge_current = False
        self.shannon_entropy = False
        
        