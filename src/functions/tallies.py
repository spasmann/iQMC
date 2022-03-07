# -*- coding: utf-8 -*-
import numpy as np

class Tallies:
    def __init__(self, init_data):
        
        self.avg_scalar_flux = init_data.avg_scalar_flux
        self.edge_scalar_flux = init_data.edge_scalar_flux
        self.avg_angular_flux = init_data.avg_angular_flux
        self.avg_current = init_data.avg_current
        self.edge_current = init_data.edge_current
        self.shannon_entropy = init_data.shannon_entropy
        
        self.Nr = init_data.Nx
        self.G = init_data.G
        
        self.phi_avg = np.zeros((self.Nr, self.G))
        self.phi_avg_old = np.zeros((self.Nr, self.G))
        self.dphi = np.zeros((self.Nr, self.G))
        self.phi_edge = np.zeros((self.Nr+1, self.G))
        self.J_avg = np.zeros((self.Nr, self.G))
        self.J_edge = np.zeros((self.Nr+1, self.G))
        self.delta_flux = 1.0
        
    def Tally(self, particle, material, mesh):
        if (self.avg_scalar_flux):
            self.AvgScalarFlux(particle, material, mesh)
            
        #if (self.avg_angular_flux):
        #    self.AvgAngularFlux()
        """
        if (self.edge_scalar_flux):
            self.EdgeScalarFlux()
        if (self.avg_current):
            self.AvgCurrent()
        if (self.edge_current):
            self.EdgeCurrent()
        """
        
    def AvgScalarFlux(self, particle, material, geometry):
        zone = particle.zone
        G = material.G
        weight = particle.weight
        ds = particle.ds
        sigt = material.sigt[zone]
        dV = geometry.CellVolume(zone)
        self.phi_avg[zone, :] += weight*(1-np.exp(-(ds*sigt)))/(sigt*dV)
    #def AvgAngularFlux(self, particle, material, geometry):
        

    def DeltaFlux(self):
        self.delta_flux = np.linalg.norm(self.phi_avg - self.phi_avg_old, np.inf)
        
    def ResetPhiAvg(self):
        self.phi_avg = np.zeros((self.Nr, self.G))
        return self.phi_avg
        