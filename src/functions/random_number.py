# -*- coding: utf-8 -*-

class RNG:
    def __init__(self, generator):
        self.generator = generator
    def newNumber(a,b):
        if (self.generator == "random"):
            return np.random.random()
    def newArray(a, b, N):
        if (self.generator == "random"):
            return np.random.uniform(a,b,N)
    def newMatrix(a, b, N1, N2):
        if (self.generator == "random"):
            return np.random.uniform(a,b,[N1,N2])