# -*- coding: utf-8 -*-

class SourceIteration:
    def __init__(init_data):
        self.init_data = init_data
        self.mesh = self.init_data.mesh
        self.material = self.init_data.material
        self.itt = 0
        self.max_iter = 50
        self.tol = 1e-5
        self.delta_flux = 1.0
        self.norm_hist = []
        self.tallies = Tallies(self.init_data.Nr, self.init_data.G)
        self.phi_avg = self.tallies.phi_avg
        self.sweep = Sweep(self.init_data,
                           self.mesh,
                           self.material, 
                           self.phi_avg,
                           self.tallies)
    def Run(self):
        while (self.itt<self.max_iter) and (self.delta_flux > self.tol):
            self.sweep.Run()
            self.tallies.DeltaFlux() 
            self.itt += 1
            self.norm_hist.append(self.tallies.delta_flux)
            self.tallies.phi_avg_old = self.tallies.phi_avg
            print("**********************")
            print("Iteration:", self.itt, "change: ",self.delta_flux)
            
    