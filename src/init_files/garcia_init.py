# -*- coding: utf-8 -*-
from mesh import Mesh

class GarciaInit:
    def __init__(self, N=100, Nx=10, generator="random"):
        self.N = N
        self.Nx = Nx
        self.generator = generator
        self.totalDim = 2
        self.RB = 5
        self.LB = 0
        self.right = False
        self.left = True
        self.material_code = "test_data"
        self.geometry = "slab"
        self.mesh = Mesh(self.Nx, self.RB)
        self.material = Material(self.material_code, self.geometry)
        
        