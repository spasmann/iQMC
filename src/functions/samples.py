# -*- coding: utf-8 -*-

import numpy as np
import math
from scipy.stats.qmc import Sobol, Halton
from src.functions.particle import Particle

class Samples:
    """
    Class for generating list of Particles given various initial conditions.
    
    This will need to become more complex to return samples for different
    sources and sampling distributions.
    """
    def __init__(self, init_data, geometry, mesh):
        self.generator = init_data.generator
        self.geometry = geometry
        self.mesh = mesh
        self.N = init_data.N
        self.Nx = init_data.Nx
        self.totalDim = init_data.totalDim
        self.RB = init_data.RB
        self.LB = init_data.LB
        self.left = init_data.left
        if (self.left):
            self.phi_left = init_data.phi_left
        self.right = init_data.right
        if (self.right):
            self.phi_right = init_data.phi_right

    def GenerateParticles(self, q):
        self.q = q
        self.counter = 0
        self.GetRnMatrix()
        self.particles = []
        self.VolumetricParticles()
        self.counter += 2 # used  to index the random number matrix 
        if (self.left):
            self.LeftBoundaryParticles()
            self.counter += 1
        if (self.right):
            self.GetRnMatrix()
            self.RightBoundaryParticles()
            self.counter += 1

    def VolumetricParticles(self):
        for i in range(self.N):
            randPos = self.rng[i,self.counter]
            randMu = self.rng[i,self.counter+1]
            pos = self.GetPos(randPos)
            mu = self.GetDir(randMu)
            zone = self.mesh.GetZone(pos)
            weight = self.VolumetricWeight(zone)
            particle = Particle(pos, mu, weight)
            self.particles.append(particle)
            
    def RightBoundaryParticles(self):
        for i in range(self.N):
            randMu = self.rng[i,self.counter]
            pos = np.array((self.RB - 1e-9,))
            mu = -np.sqrt(randMu)
            weight = self.BoundaryWeight(self.phi_right)
            particle = Particle(pos, mu, weight)
            self.particles.append(particle)
            
    def LeftBoundaryParticles(self):
        for i in range(self.N):
            randMu = self.rng[i,self.counter]
            pos = np.array((self.LB + 1e-9,))
            mu = np.sqrt(randMu)
            weight = self.BoundaryWeight(self.phi_left)
            particle = Particle(pos, mu, weight)
            self.particles.append(particle)
        
    def RandomMatrix(self):
        return np.random.uniform(0,1,[self.N,self.totalDim])
    
    def SobolMatrix(self):
        sampler = Sobol(d=self.totalDim,scramble=True)
        m = round(math.log(self.N, 2))
        return sampler.random_base2(m=m)
    
    def HaltonMatrix(self):
        sampler = Halton(d=self.totalDim,scramble=True)
        return sampler.random(n=self.N)
    
    def GetRnMatrix(self):
        if (self.generator == "random"):
            self.rng = self.RandomMatrix()
        elif (self.generator == "sobol"):
            self.rng = self.SobolMatrix()
        elif (self.generator == "halton"):
            self.rng = self.HaltonMatrix()
    
    def GetPos(self, randPos):
        return ((self.RB-self.LB)*randPos)
    
    def GetDir(self, randMu):
        return (2*randMu - 1)
    
    def GetR(self,pos):
        if (pos.size > 1):
            return np.sqrt(sum(pos**2))
        else:
            return np.abs(pos)
    
    def VolumetricWeight(self, zone):
        return self.q[zone,:]*self.geometry.CellVolume(zone)*self.Nx/self.N
    
    def BoundaryWeight(self, BV):
        # BV: boundary value, i.e. phi_left or phi_right
        return BV*self.geometry.SurfaceArea()*self.Nx/self.N
 
    

        
        
    
        