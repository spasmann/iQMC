#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue Mar 15 10:39:27 2022

@author: sampasmann
"""

import pytest
import numpy as np
from src.functions.material import Material
from src.functions.material import MaterialAvail
from src.functions.mesh import Mesh
from src.init_files import *

def material_test():
    """
    This function creates an instance of all materials available in material.py
    and runs each through a series of tests.
    """
    materials = MaterialAvail()
    geometry = "slab"
    Nx = 5
    LB = 0
    RB = 5
    mesh = Mesh(LB,RB,Nx)
    for material_code in materials:
        material = Material(material_code, geometry, mesh)
        # begin tests
        are_numpy_arrays(material)
        are_correct_size(material, Nx)
        #sum_to_sigt(material)
        matrices_in_correct_order()
    return  
    
def are_numpy_arrays(material):
    assert type(material.sigt).__module__ == np.__name__
    assert type(material.siga).__module__ == np.__name__
    assert type(material.sigs).__module__ == np.__name__
    
def are_correct_size(material, Nx):
    assert (material.sigt.shape == (Nx, material.G))
    #assert (material.siga.shape == (Nx, material.G))
    assert (material.sigs.shape == (Nx, material.G, material.G))

def sum_to_sigt(material):
    siga = material.siga
    sigs = np.sum(material.sigs,axis=0) # if multigroup the groups need to be summed
    sigt = material.sigt # need to change how this is done later
    assert (siga + sigs == sigt).all()
    return

def matrices_in_correct_order(material):
    

if __name__ == "__main__":
    material_test()
    
    