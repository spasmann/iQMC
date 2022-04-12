# -*- coding: utf-8 -*-

import numpy as np
import math
from scipy.stats.qmc import Sobol, Halton
from src.functions.particle import Particle
from src.functions.moment_matching import shift_samples

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
        self.G = init_data.G
        self.N = init_data.N
        self.Nx = init_data.Nx
        self.totalDim = init_data.totalDim
        self.RB = init_data.RB
        self.LB = init_data.LB
        self.left = init_data.left
        self.right = init_data.right
        self.moment_match = init_data.moment_match
        if (self.left):
            self.phi_left = init_data.phi_left
        if (self.right):
            self.phi_right = init_data.phi_right

    def GenerateParticles(self, q):
        self.q = q
        self.counter = 0
        self.GetRnMatrix()
        self.particles = []
        if (self.left):
            self.LeftBoundaryParticles()
            self.counter += 1
        if (self.right):
            self.RightBoundaryParticles()
            self.counter += 1
        self.VolumetricParticles()
        self.counter += 2 # used  to index the random number matrix 
        
        if (self.moment_match):
            self.moment_matching()

    def VolumetricParticles(self):
        for i in range(self.N):
            randPos = self.rng[i,self.counter]
            randMu = self.rng[i,self.counter+1]
            pos = self.GetPos(randPos) 
            mu = self.GetDir(randMu) 
            zone = self.mesh.GetZone(pos, mu)
            weight = self.VolumetricWeight(zone)
            particle = Particle(pos, mu, weight)
            self.particles.append(particle)
            
    def RightBoundaryParticles(self):
        for i in range(self.N):
            randMu = self.rng[i,self.counter]
            pos = np.array((self.RB - 1e-9,))
            mu = -np.sqrt(randMu) + 1e-9
            weight = self.BoundaryWeight(self.phi_right)
            particle = Particle(pos, mu, weight)
            self.particles.append(particle)
            
    def LeftBoundaryParticles(self):
        for i in range(self.N):
            randMu = self.rng[i,self.counter]
            pos = np.array((self.LB + 1e-9,))
            mu = np.sqrt(randMu) + 1e-9
            weight = self.BoundaryWeight(self.phi_left)
            particle = Particle(pos, mu, weight)
            self.particles.append(particle)
        
    def RandomMatrix(self):
        np.random.seed(12345)
        return np.random.uniform(0,1,[self.N,self.totalDim])
    
    def SobolMatrix(self):
        sampler = Sobol(d=self.totalDim,scramble=True)
        m = round(math.log(self.N, 2))
        return sampler.random_base2(m=m)
    
    def HaltonMatrix(self):
        sampler = Halton(d=self.totalDim,scramble=False)
        return sampler.random(n=self.N)
    
    def GetRnMatrix(self):
        if (self.generator == "random"):
            self.rng = self.RandomMatrix()
        elif (self.generator == "sobol"):
            self.rng = self.SobolMatrix()
        elif (self.generator == "halton"):
            self.rng = self.HaltonMatrix()
    
    def GetPos(self, randPos):
        return ((self.RB-self.LB)*randPos + self.LB)
    
    def GetDir(self, randMu):
        return (2*randMu - 1)
    
    def GetR(self,pos):
        if (pos.size > 1):
            return np.sqrt(sum(pos**2))
        else:
            return np.abs(pos)
    
    def VolumetricWeight(self, zone):
        weight = self.q[zone,:]*self.geometry.CellVolume(zone)/self.N*self.Nx
        return weight
    
    def BoundaryWeight(self, BV):
        # BV: boundary value, i.e. phi_left or phi_right
        weight = BV/self.N*self.geometry.SurfaceArea()
        return weight
 
    def moment_matching(self):
        ## Currently only shifting volumetric particles
        ## could shift boundary particle angle in the future
        x = np.zeros(self.N)
        mu = np.zeros(self.N)
        # we only want to shift the volumetric particles not the boundary
        start = 0
        end = self.N
        if (self.left):
            start += self.N
            end += self.N
        if (self.right):
            start += self.N
            end += self.N
        # take angle and position from voluemtric particles into new arrays
        count = 0
        for i in range(start,end):
            x[count] = self.particles[i].pos
            mu[count] = self.particles[i].dir
            count += 1
        # shift new arrays
        shifted_x = shift_samples(self.LB, self.RB, x)
        shifted_mu = shift_samples(-1.0, 1.0, mu)
        # put shifted values back into particles
        count = 0
        for j in range(start, end):
            self.particles[j].pos = shifted_x[count]
            self.particles[j].dir = shifted_mu[count]
            count += 1

        
        
    
        