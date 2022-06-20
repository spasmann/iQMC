#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue Jun  7 14:22:43 2022

@author: sampasmann
"""
import numpy as np

def GetSource(phi_avg, qmc_data):
    """
    Parameters
    ----------
    phi_avg : scalar flux matrix of size (Nx,G)
    qmc_data : an object from one of the init_files

    Returns
    -------
    q : source term, calculated per spatial cell, for source iteration

    """
    """
    try:
        assert(phi_avg.size == (Nx,G))
    except:
        print(phi_avg.size)
        Nx = qmc_data.Nx
        G = qmc_data.G
        phi_avg = np.reshape(phi_avg,(Nx,G))
    """
    material = qmc_data.material
    source   = qmc_data.source
    # calculate source for every cell individually
    q = np.zeros((material.Nx, material.G))
    for cell in range(material.Nx):
        q[cell,:] = (np.dot(phi_avg[cell,:],np.transpose(material.sigs[cell,:,:]))+source[cell,:]) 
    return q


#def GetSource(phi_avg_s, phi_avg_f, qmc_data):
    """
    GetSource(self, phi_avg_s, phi_avg_f)
    --------------------------------------
    Calculate source term for Power Iteration eigenvalue problem
    for every cell individually (loop)
    
    Parameters
    ----------
    phi_avg_s : scalar flux for inner scatter source iteration
    phi_avg_f : scalar flux for outter fission source iteration

    Returns
    -------
    q : source term



   
    q = np.empty((self.material.Nx, self.material.G), np.float64)
    for cell in range(self.material.Nx):
        q[cell,:] = (np.dot(phi_avg_s[cell,:],self.material.sigs[cell,:,:]) 
                    + np.dot(phi_avg_f[cell,:]*self.material.sigf[cell,:]*self.material.nu[cell,:],self.material.chi[cell,:])/self.k)
        #assert (np.dot(phi_avg_s[cell,:],self.material.sigs[cell,:,:])[0] >=  phi_avg_s[cell,:][0])
        #assert (np.sum(np.dot(phi_avg_s[cell,:],self.material.sigs[cell,:,:])) == np.sum(phi_avg_s[cell,:]))
    return q

    """