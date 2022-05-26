# -*- coding: utf-8 -*-
import numpy as np

class Mesh:
    def __init__(self, LB, RB, Nx):
        self.Nx = Nx
        self.dx = (RB - LB)/Nx
        self.lowR = np.linspace(LB,(RB-self.dx),Nx)
        self.highR = np.linspace(LB+self.dx,RB,Nx)
        self.midpoints = np.linspace(LB+self.dx, RB-self.dx, Nx)
    def GetZone(self, r, mu):
        if (mu > 0):    
            return np.argmax((r >= self.lowR)*(r < self.highR))
        else:
            return np.argmax((r > self.lowR)*(r <= self.highR))

    