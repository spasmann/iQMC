# -*- coding: utf-8 -*-

class Sweep:
    def __init__(self, init_data, material, phi_avg, tallies):
        self.init_data = init_data
        self.phi_avg = phi_avg
        self.material = material
        self.tallies = tallies
        self.totalDim = self.GetDim()
        self.rng = RNG(self.totalDim, self.init_data)
        self.q = GetSource()
        self.particles = Samples.GenerateParticles(self.init_data, 
                                                   self.material, 
                                                   self.q)
                
    def Run(self):
       for particle in self.particles:
           while (particle.alive)
               tallies.Tally(self.particle, self.material, self.mesh)
               particle.Move()
    
    def GetSource(self):
        if (material.G > 1):
            return self.phi_avg*np.transpose(material.sigs) + init_data.source
        else:
            return self.phi_avg*material.sigs + init_data.source
        
        
    def GetDim(self):
       
        