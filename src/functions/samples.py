# -*- coding: utf-8 -*-

import numpy as np
from particle import Particle

class Samples:
    """
    Class for generating list of Particles given various initial conditions.
    
            This will need to become more complex to return samples for different
            sources and sampling distributions.
    """
    def __init__(self, init_data):
        self.generator = init_data.generator
        self.N = init_data.N
        self.Nr = init_data.Nr
        self.totalDim = init_data.totalDim
        self.RB = init_data.RB
        self.LB = init_data.LB
        self.GetRnMatrix()
    
    def GenerateParticles(self, init_data, geometry, mesh, Q):
        self.particles = []
        for i in range(self.N):
            randPos = self.rng[0]
            randMu = self.rng[1]
            pos = self.GetPos(randPos)
            mu = self.GetDir(randMu)
            R = self.GetR(pos)
            zone = mesh.GetZone(self.)
            q = Q[zone,:]
            weight = self.GetWeight(q, zone)
            particle = Particle(pos, mu, weight)
            self.particles.append(particle)
            
    def RandomMatrix(self):
        return np.random.uniform(0,1,[self.N,self.totalDim])
    
    def GetRnMatrix(self):
        if (self.generator == "random"):
            self.rng = self.RandomMatrix()
    
    def GetPos(self, randPos):
        return ((self.RB-self.LB)*randPos)
    
    def GetDir(self, randMu):
        return (2*randMu - 1)
    
    def GetR(self,R):
        return np.sqrt(sum(pos**2))
    
    def GetWeight(self, q, zone):
        return q/self.N*geometry.CellVolume(zone)*self.Nr
    
    #def SobolMatrix():
    

        
        
    
        