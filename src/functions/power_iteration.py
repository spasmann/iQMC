# -*- coding: utf-8 -*-

from src.functions.tallies import Tallies
from src.functions.source_iteration import SourceIteration
from src.functions.sweep import Sweep
from src.functions.save_data import SaveData
import numpy as np
import matplotlib.pyplot as plt

from numba import njit

class PowerIteration:
    def __init__(self, init_data, k = 1.0):
        self.init_data = init_data
        self.source = self.init_data.source
        self.mesh = init_data.mesh
        self.material = self.init_data.material
        self.itt = 0
        self.max_SI_iter = 50
        self.max_PI_iter = 50
        self.k_tol = 1e-5
        self.phi_tol = 1e-5
        #self.norm_hist = np.empty((0,self.init_data.G))
        self.tallies = Tallies(self.init_data)
        self.sweep = Sweep(self.init_data, self.mesh, self.material)
        self.error = np.empty((0,1))
        self.k = k
        self.k_old = k+1e-4
        self.DeltaK() # difference between iterations of k (k and k_old)
        self.phi_f = np.random.random(size=(self.init_data.Nx,self.init_data.G))
    def Run(self):
        print("--------- Power Iteration ---------")
        print("Material: ", self.init_data.material_code)
        print("Random Number Generator: ", self.init_data.generator)
        print("Number of Particles per Iteration: ", self.init_data.N)
        print("Number of Spatial Cells: ", self.init_data.Nx)
        print("Initial K: ", self.k)
        # iterate over k effective
        while (self.itt<self.max_PI_iter) and (self.dk > self.k_tol):
            self.phi_f[:] = self.tallies.phi_avg[:]
            count = 0
            # iterate over scattering source
            while (count < self.max_SI_iter) and (self.tallies.delta_flux > self.phi_tol):
                self.tallies.phi_avg_old[:] = self.tallies.phi_avg[:] # shallow copy
                self.q = self.GetSource(self.tallies.phi_avg, self.phi_f)
                self.sweep.Run(self.tallies, self.q)
                self.tallies.DeltaFlux() 
                print("        Source Iteration:", count+1, "dPhi: ",self.tallies.delta_flux)
                count += 1
            self.tallies.delta_flux = 1.0
            self.UpdateK()
            self.DeltaK()
            self.itt += 1
            #self.norm_hist = np.append(self.norm_hist, self.SI.tallies.delta_flux)
            print("**********************")
            print("Power Iteration:", self.itt, "dk: ",self.dk)
            print("k: ", self.k)
            """
            if (self.init_data.true_flux.any()):
                relError = np.abs(self.tallies.phi_avg - self.init_data.true_flux)/self.init_data.true_flux
                infNorm = np.linalg.norm(relError, np.inf)
                self.error = np.append(self.error, infNorm)
            """
        if (self.init_data.save_data):
            SaveData(self.init_data, self)        


    def GetSource(self, phi_avg_s, phi_avg_f):
        """
        Calculate source term for every cell individually (loop)
        
        Parameters
        ----------
        phi_avg_s : scalar flux for inner scatter source iteration
        phi_avg_f : scalar flux for outter fission source iteration

        Returns
        -------
        q : source term

        """
        q = np.empty((self.material.Nx, self.material.G), np.float64)
        for cell in range(self.material.Nx):
            q[cell,:] = (np.dot(phi_avg_s[cell,:],self.material.sigs[cell,:,:]) 
                        + np.dot(phi_avg_f[cell,:]*self.material.sigf[cell,:]*self.material.nu[cell,:],self.material.chi[cell,:])/self.k)
            #assert (np.dot(phi_avg_s[cell,:],self.material.sigs[cell,:,:])[0] >=  phi_avg_s[cell,:][0])
            #assert (np.sum(np.dot(phi_avg_s[cell,:],self.material.sigs[cell,:,:])) == np.sum(phi_avg_s[cell,:]))
        return q
    
    
    def UpdateK(self):
        self.k_old = self.k
        self.k = self.k*(np.sum(self.material.nu*self.material.sigf*self.tallies.phi_avg)
                        /np.sum(self.material.nu*self.material.sigf*self.phi_f))

    def DeltaK(self):
        self.dk = abs(self.k-self.k_old)
        
        
        
        
    """
    def GetSource(self, phi_avg):
        # depending on the dimensionality of they problem, the source may be calculated
        # for all cells in one line, otherwise, it is calculated individually
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
    """
        
        
        
    