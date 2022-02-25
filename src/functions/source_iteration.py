# -*- coding: utf-8 -*-

from src.functions.tallies import Tallies
from src.functions.sweep import Sweep

class SourceIteration:
    def __init__(self, init_data):
        self.init_data = init_data
        self.mesh = init_data.mesh
        self.material = self.init_data.material
        self.itt = 0
        self.max_iter = 50
        self.tol = 1e-6
        self.norm_hist = []
        self.tallies = Tallies(self.init_data.Nx, self.init_data.material.G)
        self.sweep = Sweep(self.init_data,
                           self.mesh,
                           self.material)
    def Run(self):
        while (self.itt<self.max_iter) and (self.tallies.delta_flux > self.tol):
            self.tallies.phi_avg_old[:] = self.tallies.phi_avg[:] # shallow copy
            self.sweep.Run(self.tallies)
            self.tallies.DeltaFlux() 
            self.itt += 1
            self.norm_hist.append(self.tallies.delta_flux)
            print("**********************")
            print("Iteration:", self.itt, "change: ",self.tallies.delta_flux)
            
    