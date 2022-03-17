# -*- coding: utf-8 -*-
import numpy as np

class Mesh:
    def __init__(self, Nx, R):
        self.Nx = Nx
        self.dx = R[-1]/Nx
        self.lowR = np.linspace(0,(R[-1]-self.dx),Nx)
        self.highR = np.linspace(self.dx,(R[-1]),Nx)
        self.midpoints = np.linspace(self.dx, R[-1]-self.dx, Nx)
    def GetZone(self, r):
        return (np.fabs(r) >= self.lowR)*(np.fabs(r) < self.highR)
    