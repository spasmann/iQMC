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
    material = qmc_data.material
    source   = qmc_data.source
    # calculate source for every cell individually
    q = np.zeros((material.Nx, material.G))
    for cell in range(material.Nx):
        q[cell,:] = (np.dot(phi_avg[cell,:],np.transpose(material.sigs[cell,:,:]))+source[cell,:]) 
    return q