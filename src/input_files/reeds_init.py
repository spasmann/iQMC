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
    def __init__(self, N=2**12, Nx=180, generator="sobol", LB=-8.0, RB=8.0):
        self.N                  = N
        self.Nx                 = Nx
        self.generator          = generator
        self.totalDim           = 2
        self.RB                 = RB
        self.LB                 = LB
        self.G                  = 1
        self.rng_seed           = 12345
        self.RQMC               = False
        self.right              = False
        self.left               = False
        self.material_code      = "reeds_data"
        self.geometry           = "slab"
        self.mode               = "fixed_source"
        self.flux               = True
        self.flux_derivative    = False
        self.source_tilt        = False
        self.save_data          = False
        self.moment_match       = False
        self.true_flux          = reeds_sol(self.Nx,LB=LB,RB=RB)
        self.mesh               = Mesh(self.LB, self.RB, self.Nx)
        self.material           = Material(self.material_code, self.geometry, self.mesh)
        self.fixed_source       = self.material.source
        self.tallies            = Tallies(self)
        
        