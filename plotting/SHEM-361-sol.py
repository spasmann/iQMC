import numpy as np
import matplotlib.pyplot as plt

# Load XS
with np.load('~/Users/sampasmann/Documents/GitHub/QMC1D/src/materials/SHEM-361.npz') as data:
    data     = np.load('SHEM-361.npz')
    SigmaT   = data['SigmaT']
    nuSigmaF = data['nuSigmaF']
    SigmaS   = data['SigmaS']
    E        = data['E']
# Add leakage to make it subcritical
SigmaT += 0.0036 # k = 0.99112

G = len(SigmaT) # number of groups
E_mid = 0.5*(E[1:]+E[:-1]) # energy midpoint
dE    = E[1:]-E[:-1] # energy width of groups

A = np.diag(SigmaT) - SigmaS - nuSigmaF # absorption
Q = np.zeros(G) # source
Q[-2] = 1.0

phi_exact = np.linalg.solve(A,Q)/dE*E_mid
    
plt.plot(E_mid,phi_exact,label='analytical')
plt.xscale('log')
plt.xlabel(r'$E$, eV')
plt.ylabel(r'$E\phi(E)$')
plt.grid()
plt.show()
