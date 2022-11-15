# -*- coding: utf-8 -*-
import numpy as np


# =============================================================================
# Tally class
# =============================================================================
class Tallies:
    def __init__(self, init_data):
        
        self.flux               = init_data.flux
        self.flux_derivative    = init_data.flux_derivative
        self.source_tilt        = init_data.source_tilt
        self.Nr                 = init_data.Nx
        self.G                  = init_data.G
        self.dtype              = np.float64
        self.delta_flux         = 1.0
        
        if (self.flux):
            self.phi_avg = np.random.uniform(size=(self.Nr,self.G))
            self.phi_avg_old = np.random.uniform(size=(self.Nr,self.G))
        if (self.source_tilt):
            self.dphi = np.zeros((self.Nr, self.G), self.dtype)
            
    def Tally(self, particle, material, geometry, mesh):
        if (self.flux):
            avg_scalar_flux(self.phi_avg, particle, material, geometry)
        if (self.source_tilt):
            avg_scalar_flux_derivative(self.dphi, particle, material, geometry, mesh)
        #if (self.source_tilt):
        #    self.phi_avg += self.dphi#*(particle.pos[0] - mesh.midpoints[particle.zone])
            
    def DeltaFlux(self):
        self.delta_flux = np.linalg.norm(self.phi_avg - self.phi_avg_old, np.inf)
        
    def ResetPhiAvg(self):
        self.phi_avg = np.zeros((self.Nr, self.G))
        self.dphi    = np.zeros((self.Nr, self.G))
        return self.phi_avg

# =============================================================================
# Tally Functions
# =============================================================================
def avg_scalar_flux(phi_avg, particle, material, geometry):
    zone    = particle.zone
    G       = material.G
    weight  = particle.weight
    ds      = particle.ds
    sigt    = material.sigt[zone,:]
    sigt    = np.reshape(sigt, (1,G))
    dV      = geometry.CellVolume(zone)
    if (sigt.all() > 1e-12):
        phi_avg[zone,:] += (weight*(1-np.exp(-(ds*sigt)))/(sigt*dV))[0,:]
    else:
        phi_avg[zone,:] += (weight*ds/dV)    
        
def avg_scalar_flux_derivative(dphi, particle, material, geometry, mesh):
    zone    = particle.zone
    mu      = particle.angles[0]
    x       = particle.pos[0]
    x_mid   = mesh.midpoints[zone]
    dx      = mesh.dx
    G       = material.G
    w       = particle.weight
    ds      = particle.ds
    sigt    = material.sigt[zone,:]
    sigt    = np.reshape(sigt, (1,G))
    dV      = geometry.CellVolume(zone)
    if (sigt.all() > 1e-12):
        dphi[zone,:] += ((mu*(w*(1-(1+ds*sigt)*np.exp(-sigt*ds))/sigt**2) 
                        + (x - x_mid)*(w*(1-np.exp(-sigt*ds))/sigt)))[0,:]
    else:
        dphi[zone,:] += (mu*w*ds**(2)/2 + w*(x - x_mid)*ds)
    

        