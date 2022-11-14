# -*- coding: utf-8 -*-

from src.functions.geometry import Geometry
from src.functions.samples import Samples

class Sweep:
    def __init__(self, init_data):
        self.init_data  = init_data
        self.Nx         = self.init_data.Nx
        self.N          = self.init_data.N
        self.mesh       = init_data.mesh
        self.material   = init_data.material
        self.totalDim   = self.init_data.totalDim 
        self.geometry   = Geometry(init_data.geometry, self.mesh)
        self.samples    = Samples(self.init_data, self.geometry, self.mesh)
    
    def Run(self, tallies, q):
        count = 0
        self.q = q
        tallies.ResetPhiAvg()
        self.samples.GenerateParticles(self.q)
        #print("New Sweep ****************************************")
        for particle in self.samples.particles:
            particle.zone = self.mesh.GetZone(particle.pos, particle.angles)
            particle.IsAlive(self.mesh, self.geometry.geometry)
            #print("New Particle ################################")
            while (particle.alive):
                #print("-------------------")
                #print("Pos: ", particle.pos)
                #print("Angles: ", particle.angles)
                #print("Zone: ", particle.zone, "|  Radius: ", particle.R, " | ", particle.pos)
                particle.ds = self.geometry.DistanceToEdge(particle)
                tallies.Tally(particle, self.material, self.geometry, self.mesh)
                sigt = self.material.sigt[particle.zone,:]
                particle.UpdateWeight(sigt)
                particle.Move(self.mesh)
                particle.IsAlive(self.mesh, self.geometry.geometry)
                if (particle.weight.all() == 0.0):
                    particle.alive = False
            count += 1
    

        
       
        