#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon Mar 28 13:43:00 2022

@author: sampasmann
"""

from src.functions.material import Material
from src.functions.mesh import Mesh
from src.functions.tallies import Tallies
from src.input_files.reeds_solution import reeds_sol

class ReedsInit:
    def __init__(self, N=2**12, Nx=128, generator="sobol", LB=-8.0, RB=8.0,
                 source_tilt=False, seed=1234, RQMC=False):
        self.N                  = N
        self.Nx                 = Nx
        self.RB                 = RB
        self.LB                 = LB
        self.generator          = generator
        self.source_tilt        = source_tilt
        self.rng_seed           = seed
        self.RQMC               = RQMC
        self.totalDim           = 2
        self.G                  = 1
        self.keff               = 0
        self.material_code      = "reeds_data"
        self.geometry           = "slab"
        self.mode               = "fixed_source"
        self.flux               = True
        self.save_data          = False
        self.moment_match       = False
        self.right              = False
        self.left               = False
        self.true_flux          = reeds_sol(self.Nx,LB=LB,RB=RB)
        self.mesh               = Mesh(self.LB, self.RB, self.Nx, self.geometry)
        self.material           = Material(self.material_code, self.geometry, self.mesh)
        self.fixed_source       = self.material.source
        self.FixedSource        = lambda x,cell: self.fixed_source[cell,:]
        self.tallies            = Tallies(self)
        self.Nt                 = int(self.Nx*self.G)
        if (self.source_tilt):
            self.Nt = int(self.Nt*2)