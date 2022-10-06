# -*- coding: utf-8 -*-

import numpy as np
import matplotlib.pyplot as plt

Time = np.genfromtxt("Time.csv", delimiter=",")
#Time = np.flip(Time,axis=1)
Iterations = np.genfromtxt("Iterations.csv", delimiter=",")
#Iterations = np.flip(Iterations,axis=1)
Keff = np.genfromtxt("Keff.csv", delimiter=",")
#Keff = np.flip(Keff,axis=1)

figsize         = (6,5)
dpi             = 300
xlabel          = 'Pre-Conditioner Sweeps'
ylabel          = 'Restart Parameter'
edgecolors      = 'k'
linewidth       = 1
cmap            = 'cividis'

plt.figure(figsize=figsize, dpi=dpi)
plt.title('Run Time')
plt.pcolormesh(Time,edgecolors=edgecolors,linewidth=linewidth)
plt.xlabel(xlabel)
plt.ylabel(ylabel)
plt.colorbar(label="Time (seconds)")

plt.figure(figsize=figsize, dpi=dpi)
plt.title('Iterations')
plt.pcolormesh(Iterations,edgecolors=edgecolors,linewidth=linewidth)
plt.xlabel(xlabel)
plt.ylabel(ylabel)
plt.colorbar()

plt.figure(figsize=figsize, dpi=dpi)
plt.title('Davidsons Figure of Merit')
Error = abs(1.0 - Keff)
FoM = 1/(Time*Error**2)
plt.pcolormesh(FoM,edgecolors=edgecolors,linewidth=linewidth)
plt.xlabel(xlabel)
plt.ylabel(ylabel)
plt.colorbar(label=r"$1/(T*E^2)$")


"""
m = np.arange(12,0,-1)
n = np.arange(1,13)
test = np.empty((12,12))
for i in range(m.size):
    for j in range(n.size):
        test[i,j] = (m[i] + n[j]/100)
"""