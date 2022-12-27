# -*- coding: utf-8 -*-

from src.functions.geometry import Geometry
from src.functions.samples import Samples

class Sweep:
    def __init__(self, qmc_data):
        self.init_data  = qmc_data
        self.Nx         = qmc_data.Nx
        self.N          = qmc_data.N
        self.mesh       = qmc_data.mesh
        self.material   = qmc_data.material
        self.totalDim   = qmc_data.totalDim 
        self.geometry   = Geometry(qmc_data.geometry, self.mesh)
        self.samples    = Samples(qmc_data, self.geometry, self.mesh)
    
    def Run(self, qmc_data):
        count = 0
        # print(qmc_data.tallies.phi_avg)
        qmc_data.tallies.ResetPhiAvg()
        self.tallies = qmc_data.tallies
        # print(qmc_data.tallies.q)
        self.samples.GenerateParticles(self.tallies.q, self.tallies.qdot)
        #print("New Sweep ****************************************")
        for particle in self.samples.particles:
            # print()
            # print("Particle ", count+1, " ################################")
            # print("Pos ", particle.pos[0], "Angle ", particle.angles[0])
            while (particle.alive):
                #print("Zone: ", particle.zone, "|  Weight: ", particle.weight)
                particle.ds = self.geometry.DistanceToEdge(particle)
                self.tallies.Tally(particle, self.material, self.geometry, self.mesh)
                sigt = self.material.sigt[particle.zone,:]
                particle.UpdateWeight(sigt)
                particle.Move(self.mesh)
                particle.IsAlive(self.mesh, self.geometry.geometry)
                if (particle.weight.all() == 0.0):
                    particle.alive = False
            count += 1
    

        
       
        