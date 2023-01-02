#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue Jun  7 14:22:43 2022

@author: sampasmann
"""
import numpy as np

def scattering_source(cell, phi, material):
    source = np.dot((material.sigs[cell,:,:]), phi[cell,:])
    return source

def fission_source(cell, phi, keff, material):
    source = np.dot(material.nu[cell,:]*material.chi[cell,:,:]*material.sigf[cell,:]/keff, phi[cell,:])
    return source

def GetLinearSource(qmc_data):
    dphi_s       = qmc_data.tallies.dphi_s
    material     = qmc_data.material
    qdot         = np.zeros((material.Nx, material.G))
    
    for cell in range(material.Nx):
        qdot[cell,:] = (scattering_source(cell, dphi_s, material))     
    if (qmc_data.mode == "eigenvalue"):
        keff    = qmc_data.keff
        dphi_f  = qmc_data.tallies.dphi_f
        for cell in range(material.Nx):
            qdot[cell,:] += fission_source(cell, dphi_f, keff, material)
            
    qmc_data.tallies.qdot = qdot

def GetSource(phi_avg_s, qmc_data,  phi_avg_f=None):
    """
    Parameters
    ----------
    phi_avg_s : scalar flux for calculating scattering source
    qmc_data : qmc data structure
    phi_avg_f : scalar flux for calculating fission source

    Returns
    -------
    q : source in every cell, shape (Nx,G)

    """
    material        = qmc_data.material
    fixed_source    = qmc_data.fixed_source
    q               = np.zeros((material.Nx, material.G))
    if (qmc_data.source_tilt):
        GetLinearSource(qmc_data)
            
    for cell in range(material.Nx):
        q[cell,:] = (scattering_source(cell, phi_avg_s, material)  
                     + fixed_source[cell,:]) 
            
    if (phi_avg_f is not None):
        keff = qmc_data.keff
        for cell in range(material.Nx):
            q[cell,:] += fission_source(cell, phi_avg_f, keff, material) 
            
    return q
