# -*- coding: utf-8 -*-

class Geometry:
    def __init__(self,init_data):
        self.mesh = Mesh()
        self.geometry = init_data.geo
        
    def DistanceToEdge(self, particle):
        if (self.geometry == "slab"):
            return SlabEdge(particle)
        
    def SlabEdge(particle):
        if (particle.dir > 0):
            return (mesh.highR[particle.zone] - particle.R)/particle.dir
        else:
            return (mesh.lowR[particle.zone] - particle.R)/particle.dir
    
    def CellVolume(self, zone):
        if (self.geometry == "slab"):
            return (self.mesh.highR[zone] - self.mesh.lowR[zone])
        
    #def CylinderEdge():
        
    #def SphereEdge():