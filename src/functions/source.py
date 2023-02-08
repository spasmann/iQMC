#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue Jun  7 14:22:43 2022

@author: sampasmann
"""
import numpy as np


def scattering_source(cell, phi, material):
    source = np.dot((material.sigs[cell, :, :]), phi[cell, :])
    return source


def fission_source(cell, phi, keff, material):

    if (keff == 0):
        return 0
    source = np.dot(material.nu[cell, :] *
                    material.chi[cell, :, :] *
                    material.sigf[cell, :] /
                    keff, phi[cell, :])
    return source


def GetLinearSource(qmc_data):
    mode = qmc_data.mode
    dphi_s = qmc_data.tallies.dphi_s
    material = qmc_data.material
    keff = qmc_data.keff
    qdot = np.zeros((material.Nx, material.G))
    if (mode == "fixed_source"):
        dphi_f = dphi_s

    for cell in range(material.Nx):
        qdot[cell, :] = (scattering_source(cell, dphi_s, material)
                         + fission_source(cell, dphi_f, keff, material))

    qmc_data.tallies.qdot = qdot


def GetSource(phi_avg_s, qmc_data, phi_avg_f=None):
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
    mode = qmc_data.mode
    material = qmc_data.material
    keff = qmc_data.keff
    q = np.zeros((material.Nx, material.G))
    if (qmc_data.source_tilt):
        GetLinearSource(qmc_data)
    if (mode == "fixed_source"):
        phi_avg_f = phi_avg_s

    for cell in range(material.Nx):
        q[cell, :] = (scattering_source(cell, phi_avg_s, material)
                      + fission_source(cell, phi_avg_f, keff, material))

    return q
