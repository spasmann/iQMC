# -*- coding: utf-8 -*-

class Sweep:
    def __init__(self, init_data, mesh, material, phi_avg, tallies):
        self.init_data = init_data
        self.mesh = mesh
        self.material = material
        self.phi_avg = phi_avg
        self.tallies = tallies
        self.totalDim = self.GetDim() # think i want to move GetDim somewhere else
        self.geometry = Geometry(init_data)
        self.q = self.GetSource()
        self.samples = Samples(self.init_data)
        self.samples.GenerateParticles(self.init_data, 
                                       self.geometry,
                                       self.mesh,
                                       self.q)
                
    def Run(self):
       for particle in self.samples.particles:
           while (particle.alive):
               particle.zone = mesh.GetZone(particle.R)
               particle.ds = geometry.DistanceToEdge(particle)
               tallies.Tally(self.particle, 
                             self.material, 
                             self.mesh)
               sigt = self.material.sigt[particle.zone,:]
               particle.UpdateWeight(sigt)
               particle.Move()
               particle.IsAlive(self.mesh)
    

    def GetSource(self):
        if (material.G > 1):
            return self.phi_avg*np.transpose(material.sigs) + init_data.source
        else:
            return self.phi_avg*material.sigs + init_data.source
        
        
    def GetDim(self):
       
        