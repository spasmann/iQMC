# -*- coding: utf-8 -*-
"""
Want to restructure the cross section input at some point. I think a separate
file for each material/geometry, then when running the simulation the data is 
read as an input.

Sigt and Siga cross sections need to be a matrix of size (Nx, G)
Sigs needs to be (G,G)

"""
import numpy as np
import os
from src.init_files.reeds_data import reeds_data

def MaterialAvail():
    return ["first_data", "garcia_data", 12, 70, 618, "reeds_data"]

class Material:
    def __init__(self, material_code, geometry, mesh):
        self.mesh = mesh
        self.Nx = mesh.Nx
        if (material_code == "first_data"):
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
        elif (material_code == "reeds_data"):
            self.G = 1
            self.sigt, self.sigs, self.siga, self.source = reeds_data(self.Nx)
        elif (material_code == 12 or 70 or 618):
            script_dir = os.path.dirname(__file__) #<-- absolute dir the script is in
            rel_path = "../multigroup_xs/"
            abs_file_path = os.path.join(script_dir, rel_path)
            self.G = material_code
            self.D = np.genfromtxt(abs_file_path+"D_{}G_HDPE.csv".format(self.G), delimiter=",")
            self.siga = np.genfromtxt(abs_file_path+"Siga_{}G_HDPE.csv".format(self.G), delimiter=",")
            self.sigs = np.genfromtxt(abs_file_path+"Scat_{}G_HDPE.csv".format(self.G), delimiter=",")
            self.sigs = np.flip(self.sigs,1)
            self.sigt = 1/(3*self.D)
            self.sigt = np.tile(self.sigt, (self.Nx,1))  # repeat sigt so its shape (Nx, G)

            
        #    print("Type 'help(Material)' for a list of available materials")
    
