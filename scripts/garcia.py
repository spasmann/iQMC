# -*- coding: utf-8 -*-
import sys, os
sys.path.append(os.getcwd()+"/../")
from src.input_files.garcia_init import GarciaInit
from src.solvers.fixed_source.solvers import FixedSource
import matplotlib.pyplot as plt
import numpy as np
import time

from post_process.functions.functions import SN_Sweep, MOC_Sweep, garcia_angle_bins, garcia_angular_flux_sol

if __name__ == "__main__":
    # initialize problem data
    N           = 2**11
    Nx          = 24
    generator   = "sobol"
    solver      = "LGMRES"
    source_tilt = True
    s           = 1
    data        = GarciaInit(N=N, Nx=Nx, generator=generator, 
                             source_tilt=source_tilt, s=s)
    start       = time.time()
    maxit       = 50
    tol         = 1e-3
    phi         = FixedSource(data,solver=solver, maxit=maxit, tol=tol)

    stop = time.time()
    print("time: ",stop-start)
    plt.figure(dpi=300)
    plt.plot(data.mesh.midpoints,phi)
    plt.xlabel('x')
    plt.ylabel(r'$\phi$')
    plt.title('Scalar Flux')
# =============================================================================
# Plot piecewise source
# =============================================================================
    plt.figure(figsize=(6,4),dpi=300)
    plt.title('Garcia Source')
    q    = data.tallies.q
    mesh = data.mesh.edges
    mid  = data.mesh.midpoints
    x    = np.linspace(data.LB, data.RB, num=1000)
    n    = len(x)
    conditions = [(mesh[i] <= x) & (x <= mesh[i+1]) for i in range(Nx)]
    y1 = np.piecewise(x, conditions, q)
    plt.plot(x,y1, label=r'$a_j$')
    if (source_tilt):
        qdot        = data.tallies.qdot
        y2          = np.zeros_like(x)
        for i in range(n):
            zone = data.mesh.GetZone([x[i],0,0], [0,0,0])
            x_mid = data.mesh.midpoints[zone]
            y2[i] = q[zone] + qdot[zone]*(x[i] - x_mid)
        plt.plot(x,y2,label=r'$a_j + b_j(x)$')
    for i in range(len(mesh)):
        plt.axvline(mesh[i],linestyle='--',color='black')
    plt.legend()
    plt.tight_layout()

# =============================================================================
# Angular Flux 
# =============================================================================
    sol                 = garcia_angular_flux_sol(s)
    angles              = garcia_angle_bins()
    Na2                 = angles.size
    Na                  = int(Na2 / 2)
    sol_left            = sol[:Na]
    sol_right           = sol[Na:Na2]
    out                 = SN_Sweep(angles, data)
    out_left            = out[:Na]
    out_right           = out[Na:Na2]
    
    diff = np.linalg.norm((sol-out)/sol.sum())
    print("Diff: ", diff)

# =============================================================================
# Plot angular flux solutions
# =============================================================================

plt.subplots(nrows=1,ncols=2,dpi=300,figsize=(8,4))
plt.suptitle('Angular Flux Error')

plt.subplot(121)
plt.plot(angles[:Na], sol_left, 'o--', label='Sol')
plt.plot(angles[:Na], out_left, '^-', label='iQMC')
plt.xlabel(r'$\mu$')
plt.ylabel(r'$\psi(0,\mu)$')
plt.grid()
plt.title('Left Boundary')

plt.subplot(122)
plt.plot(angles[Na:Na2], sol_right, 'o--', label='Sol')
plt.plot(angles[Na:Na2], out_right, '^-', label='iQMC')
plt.xlabel(r'$\mu$')
plt.ylabel(r'$\psi(L_x,\mu)$')
plt.grid()
plt.title('Right Boundary')

plt.legend()
plt.tight_layout()
