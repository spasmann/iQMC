# -*- coding: utf-8 -*-

class Particle:
    def __init__(self, pos, dir, weight):
        self.pos = pos
        self.R = GetRadius(self.pos)
        self.dir = dir
        self.weight = weight
        
        self.distance_to_collision = 0.0
        self.distance_to_edge = 0.0
        self.alive = True
        self.group = 1
        self.zone = Mesh.GetZone(self.R)
        self.ds = 0.0
        
    def GetRadius(self):
        return np.sqrt(sum(pos**2))
        
    def Move(self, distance):
        self.pos += distance
        self.zone = Mesh.GetZone(self.R)
        
    def UpdateWeight(self):
        sigt = material.sigt
        self.weight *= np.exp(-self.ds*sigt[self.zone,:])