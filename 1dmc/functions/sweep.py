# -*- coding: utf-8 -*-

class Sweep:
    def __init__(self, init_data, material, phi_avg, tallies):
        self.init_data = init_data
        self.phi_avg = phi_avg
        self.material = material
        self.tallies = tallies
        self.totalDim = self.GetDim()
        self.rng = RNG(self.totalDim, self.init_data)
        
        if (material.G > 1):
            self.q = self.phi_avg*np.transpose(material.sigs) + init_data.source
        else:
            self.q = self.phi_avg*np.transpose(material.sigs) + init_data.source
            
        self.particles = Samples.GenerateParticles(init_data, material, self.q)
                
    def Run(self):
       for particle in self.particles:
           tallies.Tally(particle)

          
    def GetDim(self):
        
       
       
        