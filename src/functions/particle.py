# -*- coding: utf-8 -*-
import numpy as np
from mesh import Mesh

class Particle:
    def __init__(self, pos, dir, weight):
        self.pos = pos
        self.R = self.GetRadius(self.pos)
        self.dir = dir
        self.weight = weight
        
        self.distance_to_edge = 0.0
        self.alive = True
        self.group = 1
        self.ds = 0.0
        
    def GetRadius(self, pos):
        if (type(pos) == float):
            return pos
        else:
            return np.sqrt(sum(pos**2))

    def Move(self):
        self.pos += self.ds
        
    def UpdateWeight(self, sigt):
        self.weight *= np.exp(-self.ds*sigt[self.zone,:])
        
    def UpdateZone(self):
        self.zone = Mesh.GetZone(self.R)
    
    def GetZone(self, mesh):
        self.zone = mesh.GetZone(self.R)
        
    def IsAlive(self, mesh):
        if (self.R >= mesh.highR[mesh.Nr]):
            self.alive = False
