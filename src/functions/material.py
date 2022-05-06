# -*- coding: utf-8 -*-
"""

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
            self.media = 1
            
        elif (material_code == "reeds_data"):
            from src.materials.reeds_data import reeds_data
            self.sigt, self.sigs, self.siga, self.source, self.G = reeds_data(self.Nx)
            self.media = 9
            
        elif (material_code == "URRb_H2Oa5_2_0_SL_data"):
            from src.materials.URRb_H2Oa5_2_0_SL_data import URRb_H2Oa5_2_0_SL_data
            self.sigt, self.sigs, self.sigf, self.siga, self.chi, self.nu, self.G = URRb_H2Oa5_2_0_SL_data(self.Nx)
            self.media = 2
        
        elif (material_code == "PUa_1_0_SL_data"):
            from src.materials.PUa_1_0_SL_data import PUa_1_0_SL_data
            self.sigt, self.sigs, self.siga, self.sigf, self.chi, self.nu, self.G = PUa_1_0_SL_data(self.Nx)
            self.media = 1
            
        elif (material_code == "Ua_1_0_SL_data"):
            from src.materials.Ua_1_0_SL_data import Ua_1_0_SL_data
            self.sigt, self.sigs, self.siga, self.sigf, self.chi, self.nu, self.G = Ua_1_0_SL_data(self.Nx)
            self.media = 1
        
        elif (material_code == "URRa_2_0_SL"):
            from src.materials.URRa_2_0_SL_data import URRa_2_0_SL_data
            self.sigt, self.sigs, self.sigf, self.siga, self.chi, self.nu, self.G = URRa_2_0_SL_data(self.Nx)
            self.media = 1
            
        elif (material_code == "PUa_H2O_1_0_SL"):
            from src.materials.PUa_H2O_1_0_SL_data import PUa_H2O_1_0_SL_data
            self.sigt, self.sigs, self.sigf, self.siga, self.chi, self.nu, self.G = PUa_H2O_1_0_SL_data(self.Nx)
            self.media = 2
            
        elif (material_code == 12 or 70 or 618):
            from src.materials.hdpe_data import hdpe_data
            self.sigt, self.sigs, self.siga, self.G = hdpe_data(material_code,self.Nx)
            self.media = 1

    
