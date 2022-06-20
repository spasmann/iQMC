# -*- coding: utf-8 -*-

from src.functions.geometry import Geometry
from src.functions.samples import Samples
import numpy as np


class Sweep:
    def __init__(self, init_data):
        self.init_data = init_data
        self.Nx = self.init_data.Nx
        self.N = self.init_data.N
        self.mesh = init_data.mesh
        self.material = init_data.material
        self.source = self.init_data.source
        self.totalDim = self.init_data.totalDim 
        self.geometry = Geometry(init_data.geometry, self.mesh)
        self.samples = Samples(self.init_data, self.geometry, self.mesh)
    
    def Run(self, tallies, q):
        count = 0
        self.q = q
        tallies.ResetPhiAvg()
        self.samples.GenerateParticles(self.q)
        for particle in self.samples.particles:
            particle.zone = self.mesh.GetZone(particle.R, particle.dir)
            particle.IsAlive(self.mesh)
            while (particle.alive):
                particle.ds = self.geometry.DistanceToEdge(particle)
                tallies.Tally(particle, self.material, self.geometry)
                sigt = self.material.sigt[particle.zone,:]
                particle.UpdateWeight(sigt)
                particle.Move()
                particle.IsAlive(self.mesh)
                particle.zone = self.mesh.GetZone(particle.R, particle.dir)
            count += 1
    

        
       
        