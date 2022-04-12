# -*- coding: utf-8 -*-

from src.functions.tallies import Tallies
from src.functions.sweep import Sweep
from src.functions.save_data import SaveData
import numpy as np
import matplotlib.pyplot as plt

class PowerIteration:
    def __init__(self, init_data, k = 1.0):
        self.init_data = init_data
        self.source = self.init_data.source
        self.mesh = init_data.mesh
        self.material = self.init_data.material
        self.itt = 0
        self.max_iter = 50
        self.tol = 1e-5
        self.norm_hist = np.empty((0,self.init_data.G))
        self.tallies = Tallies(self.init_data)
        self.sweep = Sweep(self.init_data, self.mesh, self.material)
        self.error = np.empty((0,1))
        self.k = k
        self.k_old = (k-0.1)
        self.DeltaK() # difference between iterations of k (k and k_old)
    def Run(self):
        print("--------- Power Iteration ---------")
        print("Material: ", self.init_data.material_code)
        print("Random Number Generator: ", self.init_data.generator)
        print("Number of Particles per Iteration: ", self.init_data.N)
        print("Number of Spatial Cells: ", self.init_data.Nx)
        while (self.itt<self.max_iter) and (self.dk > self.tol):
            self.tallies.phi_avg_old[:] = self.tallies.phi_avg[:] # shallow copy
            self.q = self.GetSource(self.tallies.phi_avg)
            self.sweep.Run(self.tallies,self.q)
            self.UpdateK()
            self.DeltaK()
            self.itt += 1
            self.norm_hist = np.append(self.norm_hist, self.tallies.delta_flux)
            print("**********************")
            print("Iteration:", self.itt, r"delta k: ",self.dk)
            print("k: ", self.k)
            if (self.init_data.true_flux.any()):
                relError = np.abs(self.tallies.phi_avg - self.init_data.true_flux)/self.init_data.true_flux
                infNorm = np.linalg.norm(relError, np.inf)
                self.error = np.append(self.error, infNorm)
        
        if (self.init_data.save_data):
            SaveData(self.init_data, self)        
            
    def GetSource(self, phi_avg):
        if (self.material.G > 1):
            if (self.material.media > 1):
                q = np.zeros((self.material.Nx, self.material.G))
                for cell in range(self.material.Nx):
                    q[cell,:] = (np.dot(phi_avg[cell,:],self.material.sigs[cell,:,:]) 
                                + np.dot(phi_avg[cell,:]*self.material.sigf[cell,:]*self.material.nu[cell,:],self.material.chi[cell,:])/self.k)
                return q
            else:
                return (np.dot(phi_avg,self.material.sigs) + np.dot((phi_avg*self.material.sigf*self.material.nu).sum(axis=1),self.material.chi)/self.k)
        else:
            return (phi_avg*self.material.sigs + phi_avg*self.material.sigf*self.material.nu/self.k)
    
    def UpdateK(self):
        self.k_old = self.k
        self.k = self.k*((self.material.nu*self.material.sigf*self.tallies.phi_avg).sum()
                        /(self.material.nu*self.material.sigf*self.tallies.phi_avg_old).sum())
    def DeltaK(self):
        self.dk = abs(self.k-self.k_old)
        
        
        
    