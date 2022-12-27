# -*- coding: utf-8 -*-

import numpy as np
from src.functions.material import Material
from src.functions.mesh import Mesh
from src.functions.tallies import Tallies

class MultiGroupInit:
    def __init__(self, numGroups=12 ,N=2**12, Nx=20, generator="sobol",
                 source_tilt=False, RQMC=False):
        
        if (numGroups != 12) and (numGroups != 70) and (numGroups != 618):
            raise ValueError("Number of groups must be 12, 70, or 618.")
        
        self.N                  = N
        self.Nx                 = Nx
        self.generator          = generator
        self.source_tilt        = source_tilt
        self.RQMC               = RQMC
        self.totalDim           = 4
        self.LB                 = 0
        self.RB                 = 50
        self.G                  = numGroups
        self.rng_seed           = 12345
        self.geometry           = "slab"
        self.mode               = "fixed_source"
        self.fixed_source       = np.ones((self.Nx,self.G))
        self.material_code      = numGroups
        self.flux               = True
        self.flux_derivative    = False
        self.save_data          = False
        self.moment_match       = False
        self.right              = True
        self.left               = True
        self.mesh               = Mesh(self.LB, self.RB, self.Nx, self.geometry)
        self.material           = Material(self.material_code, self.geometry, self.mesh)
        self.true_flux          = TrueFlux(self.material, self.fixed_source, self.Nx)
        self.phi_left           = 0.5*self.true_flux[0,:]
        self.phi_left           = np.reshape(self.phi_left, (1,self.G))
        self.phi_right          = 0.5*self.true_flux[-1,:]
        self.phi_right          = np.reshape(self.phi_right, (1,self.G))
        self.tallies            = Tallies(self)
        self.Nt                 = int(self.Nx*self.G)
        if (self.source_tilt):
            self.Nt = int(self.Nt*2)
        

def TrueFlux(material, Q, Nx):
    """
    Analytic solution for infinite medium, slab, multigroup problem.
    Returns array of (Nx, G)
    """
    true_flux = np.dot(np.linalg.inv(np.diag(material.sigt[0,:]) - material.sigs[0,:,:]),Q[0,:])
    true_flux = np.tile(true_flux, (Nx,1))
    return true_flux



        