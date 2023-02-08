# -*- coding: utf-8 -*-
import numpy as np


class Particle:
    def __init__(self, pos, angles, weight, zone):
        self.pos = pos
        self.angles = angles
        self.weight = weight
        self.zone = zone
        self.R = self.GetRadius(self.pos)
        self.alive = True
        self.ds = 0.0

    def GetRadius(self, pos):
        return np.sqrt(sum(pos**2))

    def Move(self, mesh):
        x, y, z = self.pos
        mu, muSin, phi = self.angles
        ds = self.ds
        x += ds * mu
        y += ds * muSin * np.sin(phi)
        z += ds * muSin * np.cos(phi)
        self.pos = np.array((x, y, z))
        self.R = self.GetRadius(self.pos)
        self.zone = self.GetZone(mesh)

    def UpdateWeight(self, sigt):
        self.weight = self.weight * np.exp(-self.ds * sigt)

    def UpdateZone(self, mesh):
        self.zone = mesh.GetZone(self.pos, self.angles)

    def GetZone(self, mesh):
        return mesh.GetZone(self.pos, self.angles)

    def IsAlive(self, mesh, geometry):
        if (geometry == "slab"):
            if (self.angles[0] > 0):
                if (self.pos[0] >= mesh.RB):
                    self.alive = False
            else:
                if (self.pos[0] <= mesh.LB):
                    self.alive = False
        else:
            if (self.R >= mesh.RB):
                self.alive = False
