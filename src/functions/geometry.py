# -*- coding: utf-8 -*-

import numpy as np
import math

class Geometry:
    def __init__(self, geometry, mesh):
        self.geometry = geometry
        self.mesh = mesh
        
    def DistanceToEdge(self, particle):
        if (self.geometry == "slab"):
            return SlabEdge(particle, self.mesh)
        elif (self.geometry == "cylinder") or (self.geometry == "sphere"):
            return CurviLinearEdge(particle, self.mesh, self.geometry)
        
    def CellVolume(self, zone):
        highR = self.mesh.highR[zone]
        lowR  = self.mesh.lowR[zone]
        if (self.geometry == "slab"):
            return (highR - lowR)
        elif (self.geometry == "cylinder"):
            return (math.pi*(highR**2 - lowR**2))
        elif (self.geometry == "sphere"):
            return ((4.0/3.0)*math.pi*(highR**3 - lowR**3))
    
    def SurfaceArea(self):
        if (self.geometry == "slab"):
            return 0.5

def SlabEdge(particle, mesh):
    assert (particle.angles[0] != 0.0)
    x   = particle.pos[0]
    mu  = particle.angles[0]
    if (particle.angles[0] > 0.0):
        ds = (mesh.highR[particle.zone] - x)/mu + 1e-9
    else:
        ds = (mesh.lowR[particle.zone] - x)/mu + 1e-9
    return ds 

# for now, not making this part of the class so that I can use numba later
def CurviLinearEdge(particle, mesh, geometry):
    x,y,z         = particle.pos[:]
    mu,muSin,phi  = particle.angles[:]
    r             = particle.R
    zone          = particle.zone
    IB            = mesh.lowR[zone]   # inner cell boundary
    OB            = mesh.highR[zone]  # outter cell boundary
    assert (IB < r < OB)
    if (geometry == "cylinder"):
        a,k = cylinder_parameters(x,y,z,mu,muSin,phi)
    elif (geometry == "sphere"):
        a,k = sphere_parameters(x,y,z,mu,muSin,phi)
    
    # need to check both boundaries 
    for R in [IB,OB]:
        #c = (x**2 + y**2 + z**2 - R**2)
        c = (r**2 - R**2)
        if (a != 0) and (k**2 - a*c > 0):
            d1 = (-k + np.sqrt(k**2 - a*c))/a
            d2 = (-k - np.sqrt(k**2 - a*c))/a
            if (c < 0):
                if (d1 > 0):
                    distance_to_edge = d1
                else:
                    distance_to_edge = d2
            elif (d1 > 0):
                if (d1 > d2) and (d2>0):
                    distance_to_edge = d2
                else:
                    distance_to_edge = d1
    assert(distance_to_edge >= 0)
    #distance_to_edge += 1e-9
    #x           += distance_to_edge*mu 
    #y           += distance_to_edge*muSin*np.sin(phi)
    #z           += distance_to_edge*muSin*np.cos(phi)
    #R            = np.sqrt(x**2 + y**2 + z**2)
    #assert(R<=IB or R>=OB)
    return (distance_to_edge + 1e-9)


def cylinder_parameters(x,y,z,mu,muSin,phi):
    a = muSin*math.cos(phi)**2 + muSin*math.sin(phi)**2
    k = (z*muSin*math.cos(phi) + y*muSin*math.sin(phi))
    return a,k

def sphere_parameters(x,y,z,mu,muSin,phi):
    a = 1
    k = (x*mu + y*muSin*math.sin(phi) + z*muSin*math.cos(phi))
    return a,k