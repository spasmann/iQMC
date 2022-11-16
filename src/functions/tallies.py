# -*- coding: utf-8 -*-
import numpy as np

# =============================================================================
# Talies class
# =============================================================================
class Tallies:
    def __init__(self, qmc_data):
        
        self.flux               = qmc_data.flux
        self.flux_derivative    = qmc_data.flux_derivative
        self.source_tilt        = qmc_data.source_tilt
        self.Nr                 = qmc_data.Nx
        self.G                  = qmc_data.G
        self.q                  = qmc_data.source
        self.qdot               = None
        self.delta_flux         = 1.0

        if (qmc_data.mode == "eigenvalue"):
            self.phi_f       = np.random.uniform(size=(self.Nx,self.G))
        if (self.flux):
            self.phi_avg     = np.random.uniform(size=(self.Nr,self.G))
            self.phi_avg_old = np.random.uniform(size=(self.Nr,self.G))
        if (self.source_tilt):
            self.dphi_s      = np.zeros((self.Nr, self.G))
            self.qdot        = np.random.uniform(size=(self.Nr, self.G))
            if (qmc_data.mode == "eigenvalue"):
                self.dphi_f      = np.zeros((self.Nr, self.G))

# =============================================================================
# Tallies class functions
# =============================================================================

    def Tally(self, particle, material, geometry, mesh):
        if (self.flux):
            avg_scalar_flux(self.phi_avg, particle, material, geometry)
        if (self.source_tilt):
            avg_scalar_flux_derivative(self.dphi, particle, material, geometry, mesh)
            
    def DeltaFlux(self):
        self.delta_flux = np.linalg.norm(self.phi_avg - self.phi_avg_old, np.inf)
        
    def ResetPhiAvg(self):
        self.phi_avg     = np.zeros((self.Nr, self.G))
        if (self.source_tilt):
            self.dphi_s  = np.zeros((self.Nr, self.G))
            self.dphi_f  = np.zeros((self.Nr, self.G))

# =============================================================================
# Tallies
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
    

        