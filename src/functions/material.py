# -*- coding: utf-8 -*-
"""
Want to restructure the cross section input at some point. I think a separate
file for each material/geometry, then when running the simulation the data is 
read as an input.

cross sections need to be a matrix of size (Nx, G)
"""
import numpy as np

class Material:
    def __init__(self, material_code, geometry, mesh):
        self.mesh = mesh
        self.Nx = mesh.Nx
    
        if (material_code == "test_data"):
            self.G = 1
            self.sigt = np.ones((self.Nx, self.G))
            self.sigs = 0.25*self.sigt
            self.siga = self.sigt - self.sigs
        elif (material_code == "garcia_data"):
            self.G = 1
            self.sigt = np.ones((self.Nx, self.G))
            self.c = 1.0
            self.sigs = np.exp(-self.mesh.midpoints/self.c)
            self.sigs = np.reshape(self.sigs, (self.Nx, self.G))
            self.siga = self.sigt - self.sigs
        elif (material_code == 12):
            self.G = material_code
            self.D = np.genfromtxt("../src/multigroup_xs/D_{}G_HDPE.csv".format(self.G), delimiter=",")
            self.siga = np.genfromtxt("../src/multigroup_xs/Siga_{}G_HDPE.csv".format(self.G), delimiter=",")
            self.sigs = np.genfromtxt("../src/multigroup_xs/Scat_{}G_HDPE.csv".format(self.G), delimiter=",")
            self.sigs = np.flip(self.sigs,1)
            self.sigt = 1/(3*self.D)
            # repeat sigt so its shape (Nx, G)
            self.sigt = np.tile(self.sigt, (self.Nx,1))
            
        #else:
        #    print("Type 'help(Material)' for a list of available materials")
    
        