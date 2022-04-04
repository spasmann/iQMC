# -*- coding: utf-8 -*-
"""
Want to restructure the cross section input at some point. I think a separate
file for each material/geometry, then when running the simulation the data is 
read as an input.

Sigt and Siga cross sections need to be a matrix of size (Nx, G)

"""

def MaterialAvail():
    return ["first_data", "garcia_data", 12, 70, 618, "reeds_data"]

class Material:
    def __init__(self, material_code, geometry, mesh):
        self.mesh = mesh
        self.Nx = mesh.Nx
        if (material_code == "garcia_data"):
            from src.materials.garcia_data import garcia_data
            self.sigt, self.sigs, self.siga, self.G = garcia_data(self.mesh, self.Nx)
        elif (material_code == "reeds_data"):
            from src.materials.reeds_data import reeds_data
            self.sigt, self.sigs, self.siga, self.source, self.G = reeds_data(self.Nx)
        elif (material_code == 12 or 70 or 618):
            from src.materials.hdpe_data import hdpe_data
            self.sigt, self.sigs, self.siga, self.G = hdpe_data(material_code,self.Nx)

    
