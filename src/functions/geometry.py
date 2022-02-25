# -*- coding: utf-8 -*-

class Geometry:
    def __init__(self, geometry, mesh):
        #self.mesh = Mesh()
        self.geometry = geometry
        self.mesh = mesh
        
    def DistanceToEdge(self, particle):
        if (self.geometry == "slab"):
            return self.SlabEdge(particle)
        
    def SlabEdge(self, particle):
        if (particle.dir > 0):
            return (self.mesh.highR[particle.zone] - particle.R)/particle.dir + 1e-9
        else:
            return (self.mesh.lowR[particle.zone] - particle.R)/particle.dir + 1e-9
    
    def CellVolume(self, zone):
        if (self.geometry == "slab"):
            return (self.mesh.highR[zone] - self.mesh.lowR[zone])
    
    def SurfaceArea(self):
        return 0.5
        
    #def CylinderEdge():
        
    #def SphereEdge():