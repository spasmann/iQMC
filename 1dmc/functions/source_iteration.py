# -*- coding: utf-8 -*-

class SourceIteration:
    def __init__():
        self.itt = 0
        self.max_iter = 50
        self.tol = 1e-5
        self.delta_flux
        self.norm_hist = []
        self.tallies = Tallies()
        self.phi_avg = self.tallies_phi_avg
        self.sweep = Sweep(self.init_data, 
                           self.material, 
                           self.tallies)
    def Run(self):
        while (self.itt<self.max_iter) and (self.delta_flux > self.tol):
            self.sweep.Run()
            self.tallies.DeltaFlux() 
            self.itt += 1
            self.norm_hist.append(self.delta_flux)
            print("**********************")
            print("Iteration:", self.itt, "change: ",self.delta_flux)
            
    