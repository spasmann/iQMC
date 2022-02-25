# -*- coding: utf-8 -*-
import numpy as np

class Tallies:
    def __init__(self, Nr, G):
        self.avg_scalar_flux = True
        self.edge_scalar_flux = False
        self.avg_angular_flux = False
        self.avg_current = False
        self.edge_current = False
        self.shannon_entropy = False
        
        self.Nr = Nr
        self.G = G
        self.phi_avg = np.zeros((Nr, G))
        self.phi_avg_old = np.zeros((Nr, G))
        self.dphi = np.zeros((Nr, G))
        self.phi_edge = np.zeros((Nr+1, G))
        self.J_avg = np.zeros((Nr, G))
        self.J_edge = np.zeros((Nr+1, G))
        self.delta_flux = 1.0
        
    def Tally(self, particle, material, mesh):
        if (self.avg_scalar_flux):
            self.AvgScalarFlux(particle, material, mesh)
        """
        if (self.edge_scalar_flux):
            self.EdgeScalarFlux()
        if (self.avg_angular_flux):
            self.AvgAngularFlux()
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

    def DeltaFlux(self):
        self.delta_flux = np.linalg.norm(self.phi_avg - self.phi_avg_old, np.inf)
        
    def ResetPhiAvg(self):
        self.phi_avg = np.zeros((self.Nr, self.G))
        return self.phi_avg
        