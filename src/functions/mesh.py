# -*- coding: utf-8 -*-
import numpy as np

class Mesh:
    def __init__(self, LB, RB, Nx):
        self.Nx         = Nx
        self.LB         = LB
        self.RB         = RB
        self.dx         = (RB - LB)/Nx
        self.lowR       = np.linspace(LB,(RB-self.dx),Nx)
        self.highR      = np.linspace(LB+self.dx,RB,Nx)
        self.midpoints  = np.linspace(LB+(self.dx*0.5), RB-(self.dx*0.5), Nx)
        self.edges      = np.linspace(LB,RB,Nx+1)
        
    def GetZone(self, pos, angles):
        if (pos[1] == 0) and (pos[2] == 0):
            mu = angles[0]
            x  = pos[0]
            if (mu > 0):    
                return np.argmax((x >= self.lowR)*(x < self.highR)) 
            else:
                return np.argmax((x > self.lowR)*(x <= self.highR)) 
        else:
            r = np.sqrt(sum(pos**2))
            return np.argmax((r > self.lowR)*(r < self.highR)) 
        

    