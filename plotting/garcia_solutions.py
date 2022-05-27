#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed May 11 13:00:26 2022

@author: sampasmann
"""

import numpy as np
import matplotlib.pyplot as plt

left1    = [5.89664e-01,5.31120e-01,4.43280e-01,3.80306e-01,3.32964e-01,2.96090e-01,2.66563e-01,2.42390e-01,2.22235e-01,2.05174e-01,1.90546e-01]    
right1   = [6.07488e-06,6.92516e-06,9.64232e-06,1.62339e-05,4.38580e-05,1.69372e-04,5.73465e-04,1.51282e-03,3.24369e-03,5.96036e-03,9.77123e-03]

leftInf  = [8.97798e-01,8.87836e-01,8.69581e-01,8.52299e-01,8.35503e-01,8.18996e-01,8.02676e-01,7.86493e-01,7.70429e-01,7.54496e-01,7.38721e-01]
rightInf = [1.02202e-01,1.12164e-01,1.30419e-01,1.47701e-01,1.64497e-01,1.81004e-01,1.97324e-01,2.13507e-01,2.29571e-01,2.45504e-01,2.61279e-01]

angle = np.array((0.05,0.1,0.2,0.3,0.4,0.5,0.6,0.7,0.8,0.9,1.0))

plt.figure(dpi=300,figsize=(8,4))

plt.subplot(121)
plt.plot(-angle,left1,'-o',label=r'$s=1.0$')
plt.plot(-angle,leftInf, '--o', label=r'$s=\infty$')
plt.xlabel(r'$\mu$')
plt.ylabel(r'$\psi(0,\mu)$')
plt.legend()
plt.tight_layout()
plt.subplot(122)
plt.plot(angle,right1,'-s',label=r'$s=1.0$')
plt.plot(angle,rightInf,'--s',label=r'$s=\infty$')
plt.xlabel(r'$\mu$')
plt.ylabel(r'$\psi(\tau,\mu)$')
plt.legend()
plt.tight_layout()
