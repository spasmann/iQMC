# -*- coding: utf-8 -*-
import numpy as np

class Mesh:
    def __init__(self, LB, RB, Nx):
        self.Nx = Nx
        self.LB = LB
        self.RB = RB
        self.dx = (RB - LB)/Nx
        self.lowR = np.linspace(LB,(RB-self.dx),Nx)
        self.highR = np.linspace(LB+self.dx,RB,Nx)
        self.midpoints = np.linspace(LB+(self.dx*0.5), RB-(self.dx*0.5), Nx)
    def GetZone(self, r, mu):
        # this function might not work anymore because mu isnt the only angle
        if (mu > 0):    
            return np.argmax((r >= self.lowR)*(r < self.highR)) 
        else:
            return np.argmax((r > self.lowR)*(r <= self.highR)) 

    