# -*- coding: utf-8 -*-

from src.functions.geometry import Geometry
from src.functions.samples import Samples


class Sweep:
    def __init__(self, qmc_data):
        self.init_data = qmc_data
        self.Nx = qmc_data.Nx
        self.N = qmc_data.N
        self.mesh = qmc_data.mesh
        self.material = qmc_data.material
        self.totalDim = qmc_data.totalDim
        self.geometry = Geometry(qmc_data.geometry, self.mesh)
        self.samples = Samples(qmc_data, self.geometry, self.mesh)

    def Run(self, qmc_data):
        count = 0
        qmc_data.tallies.ResetPhiAvg()
        self.tallies = qmc_data.tallies
        self.samples.GenerateParticles(self.tallies.q, self.tallies.qdot)
        for particle in self.samples.particles:
            while (particle.alive):
                particle.ds = self.geometry.DistanceToEdge(particle)
                self.tallies.Tally(
                    particle,
                    self.material,
                    self.geometry,
                    self.mesh)
                sigt = self.material.sigt[particle.zone, :]
                particle.UpdateWeight(sigt)
                particle.Move(self.mesh)
                particle.IsAlive(self.mesh, self.geometry.geometry)
                if (particle.weight.all() == 0.0):
                    particle.alive = False
            count += 1
