#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue Mar 15 10:57:01 2022

@author: sampasmann
"""

import pytest
import numpy as np
from functions.mesh import Mesh

def mesh_test():
    R = np.array((np.random.random()*10,))
    Nx = int(np.random.random()*1000)
    mesh = Mesh(Nx, R)
    """
    This function creates an instance of all problems defined in the problems/
    folder and and runs tests on all created meshes
    """
    materials = ProblemsAvail()
    geometry = "slab"
    Nx = 5
    RB = 5
    mesh = Mesh(Nx, np.array((RB,)))
    for material_code in materials:
        