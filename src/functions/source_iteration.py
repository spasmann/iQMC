# -*- coding: utf-8 -*-

from src.functions.tallies import Tallies
from src.functions.sweep import Sweep
from src.functions.save_data import SaveData
import numpy as np

class SourceIteration:
    def __init__(self, init_data):
        self.init_data = init_data
        self.save_data = self.init_data.save_data 
        self.source = self.init_data.source
        self.mesh = init_data.mesh
        self.material = self.init_data.material
        self.itt = 0
        self.max_iter = 50
        self.tol = 1e-6
        self.norm_hist = np.empty((0,self.init_data.G))
        self.tallies = Tallies(self.init_data)
        self.sweep = Sweep(self.init_data,
                           self.mesh,
                           self.material)
        self.error = np.empty((0,1))
        self.fname = "default"
    def Run(self):
        print("--------- Source Iteration ---------")
        print("Material: ", self.init_data.material_code)
        print("Random Number Generator: ", self.init_data.generator)
        print("Number of Particles per Iteration: ", self.init_data.N)
        print("Number of Spatial Cells: ", self.init_data.Nx)
        while (self.itt<self.max_iter) and (self.tallies.delta_flux > self.tol):
            self.tallies.phi_avg_old[:] = self.tallies.phi_avg[:] # shallow copy
            self.q = self.GetSource(self.tallies.phi_avg)
            self.sweep.Run(self.tallies, self.q)
            self.tallies.DeltaFlux() 
            self.itt += 1
            self.norm_hist = np.append(self.norm_hist, self.tallies.delta_flux)
            print("**********************")
            print("Iteration:", self.itt, "change: ",self.tallies.delta_flux)
            if (self.init_data.true_flux.any()):
                #relError = np.abs(self.tallies.phi_avg - self.init_data.true_flux)/self.init_data.true_flux
                Err = np.abs(self.tallies.phi_avg - self.init_data.true_flux)
                infNorm = np.linalg.norm(Err, np.inf)
                self.error = np.append(self.error, infNorm)
        
        if (self.save_data):
            SaveData(self.init_data, self, self.fname)

    """
    def GetSource(self, phi_avg):
        if (self.material.G > 1):
            return (np.dot(phi_avg,np.transpose(self.material.sigs)) + self.source)
        else:
            return (phi_avg*self.material.sigs + self.source)
    """
    
    def GetSource(self, phi_avg):
        # calculate source for every cell individually
        q = np.zeros((self.material.Nx, self.material.G))
        for cell in range(self.material.Nx):
            q[cell,:] = (np.dot(phi_avg[cell,:], (self.material.sigs[:,:])+self.source[cell,:]) )
        return q
    
    
    
    