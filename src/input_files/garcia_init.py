# -*- coding: utf-8 -*-

import numpy as np
from src.functions.material import Material
from src.functions.mesh import Mesh
from src.functions.tallies import Tallies

class GarciaInit:
    def __init__(self, N=2**12, Nx=20, generator="sobol", source_tilt=False,
                 seed=12345, RQMC=False):        
        self.N                  = N
        self.Nx                 = Nx
        self.generator          = generator
        self.source_tilt        = source_tilt
        self.rng_seed           = seed
        self.RQMC               = RQMC
        self.totalDim           = 3
        self.RB                 = 5.0
        self.LB                 = 0.0
        self.G                  = 1
        self.c                  = np.inf
        self.keff               = 0
        self.right              = False
        self.left               = True
        self.phi_left           = 1.0
        self.phi_right          = 0.0
        self.material_code      = "garcia_data"
        self.geometry           = "slab"
        self.mode               = "fixed_source"
        self.flux               = True
        self.save_data          = False
        self.moment_match       = False
        self.true_flux          = np.array((False))
        self.mesh               = Mesh(self.LB, self.RB, self.Nx, self.geometry)
        self.material           = Material(self.material_code, self.geometry, self.mesh)
        self.fixed_source       = np.zeros((self.Nx,self.G))
        self.FixedSource        = lambda x, cell: self.fixed_source[cell,:]
        self.tallies            = Tallies(self)
        self.Nt                 = int(self.Nx*self.G)
        if (self.source_tilt):
            self.Nt = int(self.Nt*2)
        
        