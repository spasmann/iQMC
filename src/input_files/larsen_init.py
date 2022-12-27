# -*- codin#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon Mar 28 13:43:00 2022

@author: sampasmann
"""

from src.functions.mesh import Mesh
from src.functions.tallies import Tallies
from src.functions.material import Material
import numpy as np

class LarsenInit:
    def __init__(self, N=2**12, Nx=20, generator="sobol", LB=0.0, RB=5.0,
                 source_tilt=False):
        self.N                  = N
        self.Nx                 = Nx
        self.RB                 = RB
        self.LB                 = LB
        self.generator          = generator
        self.source_tilt        = source_tilt
        self.totalDim           = 4
        self.G                  = 1
        self.rng_seed           = 12345
        self.q0                 = 1.0
        self.q1                 = 1.0
        self.material_code      = "larsen_data"
        self.geometry           = "slab"
        self.mode               = "fixed_source"
        self.flux               = True
        self.save_data          = False
        self.moment_match       = False
        self.RQMC               = False
        self.right              = True
        self.left               = True
        self.mesh               = Mesh(self.LB, self.RB, self.Nx, self.geometry)
        self.material           = Material(self.material_code, self.geometry, self.mesh)
        self.fixed_source       = larsen_source(self.q0, self.q1, self.G, self.mesh)
        self.phi_right          = 0.5*self.fixed_source[-1,:]
        self.phi_left           = 0.5*self.fixed_source[0,:]
        self.true_flux          = larsen_sol(self.fixed_source, self.material)
        self.tallies            = Tallies(self)
        self.Nt                 = int(self.Nx*self.G)
        if (self.source_tilt):
            self.Nt = int(self.Nt*2)
    

def larsen_source(q0, q1, G, mesh):
    midpoints       = mesh.midpoints
    Nx              = mesh.Nx
    fixed_source    = q0 + q1*midpoints
    fixed_source    = np.reshape(fixed_source, (Nx,G))
    return fixed_source

def larsen_material(mesh):
    Nx   = mesh.Nx
    G    = 1
    sigt = 1.0
    sigs = 0.2
    siga = sigt - sigs
    
    sigt = np.tile(sigt, (Nx,G))
    sigs = np.tile(sigs, (Nx,G,G))
    siga = np.tile(siga, (Nx,G))
    
    return sigt, sigs, siga, G

def larsen_sol(fixed_source, material):
    return fixed_source/material.siga