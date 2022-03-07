# -*- coding: utf-8 -*-
import numpy as np

class Mesh:
    def __init__(self, Nr, R):
        self.Nx = Nr
        self.dr = R[-1]/Nr
        self.lowR = np.linspace(0,(R[-1]-self.dr),Nr)
        self.highR = np.linspace(self.dr,(R[-1]),Nr)
        self.midpoints = np.linspace(self.dr, R[-1]-self.dr, Nr)
    def GetZone(self, r):
        return (np.fabs(r) >= self.lowR)*(np.fabs(r) < self.highR)
    
    def CellVolume(self, zone):
        if (self.geometry == "slab"):
            return (self.mesh.highR[zone] - self.mesh.lowR[zone])