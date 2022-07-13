#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue Jul 12 12:49:47 2022

@author: sampasmann
"""

from scipy.stats.qmc import Sobol
import matplotlib.pyplot as plt

m = 10
sampler = Sobol(d=2, scramble=False)
samples = sampler.random_base2(m=m)
# 1280×640px for best display
my_dpi = 300
plt.figure(figsize=(1280/my_dpi, 640/my_dpi), dpi=my_dpi)
plt.plot(samples[:,0], samples[:,1],'o',markersize=2)
plt.axis('off')
plt.savefig('sobol.png', transparent=True)
#plt.title('iQMC')
