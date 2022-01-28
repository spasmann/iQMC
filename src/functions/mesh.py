# -*- coding: utf-8 -*-

class Mesh:
    def __init__(self, Nr, R):
        self.dr = R[-1]/Nr
        self.lowR = np.linspace(0,(R[-1]-dr),Nr)
        self.highR = np.linspace(dr,(R[-1]),Nr)
    def GetZone(self, r):
        return (np.fabs(r) >= lowR)*(np.fabs(r) < highR)