# -*- coding: utf-8 -*-

class Tallies:
    def __init__(self):
        self.avg_scalar_flux = True
        self.edge_scalar_flux = False
        self.avg_angular_flux = False
        self.avg_current = False
        self.edge_current = False
        self.shannon_entropy = False
        
        self.phi_avg = np.zeros(Nr, G)
        self.phi_avg_old = np.zeros(Nr, G)
        self.dphi = np.zeros(Nr, G)
        self.phi_edge = np.zeros(Nr+1, G)
        self.J_avg = np.zeros(Nr, G)
        self.J_edge = np.zeros(Nr+1, G)
        self.delta_flux = 0.0
        
    def Tally(self, particle, material, mesh):
        if (self.avg_scalar_flux):
            self.AvgScalarFlux()
        if (self.edge_scalar_flux):
            self.EdgeScalarFlux()
        if (self.avg_angular_flux):
            self.AvgAngularFlux()
        if (self.avg_current):
            self.AvgCurrent()
        if (self.edge_current):
            self.EdgeCurrent()
            
    def AvgScalarFlux(self, particle, material, mesh):
        zone = particle.zone
        G = particle.G
        weight = particle.weight
        ds = particle.PathLength(zone)
        sigt = material.sigt(zone)
        dV = mesh.dV(zone)
        self.avg_scalar_flux[zone, G] += weight*(1-np.exp(-(ds*sigt)))/(sigt*dV)

    def DeltaFlux(self):
        self.delta_flux = np.linalg.norm(self.phi_avg - self.phi_avg_old)
        