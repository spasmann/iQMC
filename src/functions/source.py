#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue Jun  7 14:22:43 2022

@author: sampasmann
"""
import numpy as np

def scattering_source(cell, phi, material):
    source = np.dot((material.sigs[cell,:,:]), phi[cell,:])
    #source = np.zeros(material.G)
    #for i in range(material.G):
    #    source[i] = np.dot(material.sigs[cell,i,:], phi[cell,:])
    return source

def fission_source(cell, phi, keff, material):
    source = np.dot(material.nu[cell,:]*material.chi[cell,:,:]*material.sigf[cell,:]/keff, phi[cell,:])
    #source = np.zeros(material.G)
    #for i in range(material.G):
    #    source[i] = np.dot(material.chi[cell,i]*material.nu[cell,:]*material.sigf[cell,:]/keff, phi[cell,:])
    return source

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
    fixed_source    = qmc_data.source
    q               = np.zeros((material.Nx, material.G))
    if (phi_avg_f is None):
        for cell in range(material.Nx):
            q[cell,:] = (scattering_source(cell, phi_avg_s, material)  
                         + fixed_source[cell,:]) 
    else:
        keff = qmc_data.keff
        for cell in range(material.Nx):
            q[cell,:] = (scattering_source(cell, phi_avg_s, material) 
                         + fission_source(cell, phi_avg_f, keff, material) 
                         + fixed_source[cell,:]) 
    return q
