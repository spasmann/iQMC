# -*- coding: utf-8 -*-

import numpy as np
import math
from scipy.stats.qmc import Sobol, Halton, LatinHypercube
from src.functions.particle import Particle
from mpi4py import MPI

class Samples:
    """
    Class for generating list of Particles given various initial conditions.
    
    This will need to become more complex to return samples for different
    sources and sampling distributions.
    """
    def __init__(self, init_data, geometry, mesh):
        self.generator  = init_data.generator
        self.RQMC       = False
        self.geometry   = geometry
        self.mesh       = mesh
        self.G          = init_data.G
        self.N          = init_data.N
        self.Nx         = init_data.Nx
        self.totalDim   = init_data.totalDim
        self.RB         = init_data.RB
        self.LB         = init_data.LB
        self.left       = init_data.left
        self.right      = init_data.right
        
        # split the total number of particles into the different sources
        #self.Nleft = 0
        #self.Nright=0
        #self.Nvolumetric = self.N
        self.moment_match = init_data.moment_match
        if (self.left):
            self.phi_left = init_data.phi_left
            #self.left = math.floor(0.125*self.N)
        if (self.right):
            self.phi_right = init_data.phi_right
            #self.right = math.floor(0.125*self.N)
            
        #self.Nvolumetric = self.N - self.Nright - self.Nleft
        # use MPI rank and nproc to determine which random numbers to use
        # each rank will generate the whole matrix, then use "start" and 
        # "stop" to grab the appropriate section
        comm        = MPI.COMM_WORLD
        rank        = comm.Get_rank()
        nproc       = comm.Get_size()
        self.rank   = comm.Get_rank()
        self.nproc  = comm.Get_size()
        self.start  = math.floor((rank/nproc)*self.N)
        self.stop   = math.floor((rank+1)/nproc*self.N) 

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
        geo = self.geometry.geometry
        for i in range(self.start,self.stop):
            randX   = self.rng[i,self.counter]
            randMu  = self.rng[i,self.counter+1]
            x       = self.GetPos(randX) 
            mu      = self.GetMu(randMu)
            if (mu == 0.0):
                mu += 0.01
            if (geo == "cylinder") or (geo == "sphere"):
                randPhi = self.rng[i,self.counter+2]
                phi     = self.GetPhi(randPhi)
                muSin   = math.sqrt(1-mu**2)
                angle   = np.array((mu, muSin, phi))
            else:
                angle   = np.array((mu, 0, 0)) # mu, muSin, phi
            pos      = np.array((x,0,0)) # x, y, z
            zone     = self.mesh.GetZone(pos, angle)
            weight   = self.VolumetricWeight(zone)
            particle = Particle(pos, angle, weight, zone)
            self.particles.append(particle)
            
    def RightBoundaryParticles(self):
        pos = np.array((self.RB - 1e-9,))
        weight = self.BoundaryWeight(self.phi_right)
        for i in range(self.start,self.stop):
            randMu = self.rng[i,self.counter]
            mu = -np.sqrt(randMu) - 1e-9
            particle = Particle(pos, mu, weight)
            self.particles.append(particle)
            
    def LeftBoundaryParticles(self):
        pos = np.array((self.LB + 1e-9,))
        weight = self.BoundaryWeight(self.phi_left)
        for i in range(self.start,self.stop):
            randMu = self.rng[i,self.counter]
            mu = np.sqrt(randMu) + 1e-9
            particle = Particle(pos, mu, weight)
            self.particles.append(particle)
        
    def RandomMatrix(self):
        np.random.seed(56789)
        return np.random.random((self.N,self.totalDim))
    
    def SobolMatrix(self):
        sampler = Sobol(d=self.totalDim,scramble=self.RQMC, seed=1234)
        m = round(math.log(self.N, 2))
        samples = sampler.random_base2(m=m)
        return samples
    
    def HaltonMatrix(self):
        sampler = Halton(d=self.totalDim,scramble=self.RQMC)
        sampler.fast_forward(1)
        return sampler.random(n=self.N)
    
    def LatinHypercube(self):
        sampler = LatinHypercube(d=self.totalDim)
        return sampler.random(n=self.N)
    
    def GetRnMatrix(self):
        if (self.generator == "random"):
            self.rng = self.RandomMatrix()
        elif (self.generator == "sobol"):
            self.rng = self.SobolMatrix()
        elif (self.generator == "halton"):
            self.rng = self.HaltonMatrix()
        elif (self.generator == "latin_hypercube"):
            self.rng = self.LatinHypercube()
    
    def GetPos(self, randPos):
        return ((self.RB-self.LB)*randPos + self.LB)
    
    def GetMu(self, randMu):
        return (2*randMu - 1)
    
    def GetPhi(self, randPhi):
        return (2*math.pi*randPhi)
    
    def GetR(self,pos):
        return np.sqrt(sum(pos**2))

    def VolumetricWeight(self, zone):
        weight = self.q[zone,:]*self.geometry.CellVolume(zone)/self.N*self.Nx
        return weight
    
    def BoundaryWeight(self, BV):
        # BV: boundary value, i.e. phi_left or phi_right 
        return (self.geometry.SurfaceArea()*BV/self.N)