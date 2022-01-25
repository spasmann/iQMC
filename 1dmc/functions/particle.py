# -*- coding: utf-8 -*-

class Particle:
    def __init__(self, pos, dir, wgt):
        self.pos = pos
        self.dir = dir
        self.weight = weight
        
        self.distance_to_collision = 0.0
        self.distance_to_edge = 0.0
        self.collide = False
        self.is_alive = True
        self.group = 1
        self.zone = Mesh.getZone(self.pos)
        
    def getDistanceToCollision(self, RNG, material):
        rn = RNG.newNumber(0,1)
        sigt = material.sigt(self.zone, self.group) # might need to change to []
        self.distance_to_collision = -np.log(1.0-rn)/sigt
        
    def getDistanceToEdge():
        
        
    def move(self, distance):
        self.pos += distance
        
    def updateWeight(self):
        self.weight *= np.exp(-ds*sigt[self.zone,:])