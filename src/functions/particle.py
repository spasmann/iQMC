# -*- coding: utf-8 -*-
import numpy as np

class Particle:
    def __init__(self, pos, dir, weight):
        self.pos = pos
        self.weight = weight
        self.R = self.GetRadius(self.pos)
        self.dir = dir
        self.alive = True
        self.ds = 0.0
        self.zone = 0
        
    def GetRadius(self, pos):
        if (pos.size > 1):
            return np.sqrt(sum(pos**2))
        else:
            return pos

    def Move(self):
        self.pos += self.ds*self.dir
        self.R = self.GetRadius(self.pos)
        
    def UpdateWeight(self, sigt):
        self.weight *= np.exp(-self.ds*sigt)
        
    def UpdateZone(self, mesh):
        self.zone = Mesh.GetZone(self.R, self.dir)
    
    def GetZone(self, mesh):
        self.zone = mesh.GetZone(self.R, self.dir)
        
    def IsAlive(self, mesh):
        if (self.dir > 0):
            if (self.R >= mesh.RB):
                self.alive = False
        else:
            if (self.R <= mesh.LB):
                self.alive = False
