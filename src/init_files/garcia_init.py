# -*- coding: utf-8 -*-

import numpy as np
from src.functions.material import Material
from src.functions.mesh import Mesh

class GarciaInit:
    def __init__(self, N=2**12, Nx=20, generator="random"):
        self.N = N
        self.Nx = Nx
        self.generator = generator
        self.totalDim = 3
        self.RB = 5
        self.LB = 0
        self.G = 1
        self.right = False
        self.left = True
        self.phi_left = 1.0
        self.material_code = "garcia_data"
        self.geometry = "slab"
        self.mesh = Mesh(self.Nx, [self.RB])
        self.material = Material(self.material_code, self.geometry, self.mesh)
        self.source = np.zeros((self.Nx,self.G))
        
        