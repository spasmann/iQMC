#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue Jun  7 13:18:01 2022

@author: sampasmann
"""
import time
import numpy as np
from mpi4py import MPI
from src.functions.save_data import SaveData
from src.solvers.fixed_source.solvers import Picard
from src.solvers.eigenvalue.maps import MatVec_data, MatVec
from scipy.sparse.linalg import gmres, lgmres, bicgstab, LinearOperator
import scipy.linalg as sp
from src.solvers.eigenvalue.maps import SI_Map

# =============================================================================
# Iteration and Residual Storage for Krylov Solvers
# =============================================================================


class gmres_counter(object):
    def __init__(self, disp=True):
        self._disp = disp
        self.iter = 0
        self.callbacks = []

    def __call__(self, rk=None):
        self.callbacks.append(rk.copy())
        self.iter += 1
        if self._disp:
            if (self.iter > 1):
                print("     Iteration:", self.iter - 1, "change: ",
                      np.linalg.norm((rk - self.callbacks[self.iter - 2])))

# =============================================================================
# Power Iteration
# =============================================================================
# TODO: Picard PI is not working


def PowerIteration(qmc_data, solver="LGMRES", max_outter_itt=10,
                   max_inner_itt=10, outter_tol=1e-5, inner_tol=1e-5,
                   report_progress=True):
    comm = MPI.COMM_WORLD
    rank = comm.Get_rank()
    # nproc       = comm.Get_size()

    itt = 0
    k = qmc_data.keff
    dk = 1.0
    phi_old = qmc_data.tallies.phi_f.copy()

    # res_hist        = []
    k_hist = []
    if (rank == 0):
        print("")
        print("    ██╗ ██████╗ ███╗   ███╗ ██████╗")
        print("    ║ ║██╔═══██╗████╗ ████║██╔════╝")
        print("    ██║██║   ██║██╔████╔██║██║     ")
        print("    ██║██║▄▄ ██║██║╚██╔╝██║██║     ")
        print("    ██║╚██████╔╝██║ ╚═╝ ██║╚██████╗")
        print("    ╚═╝ ╚══▀▀═╝ ╚═╝     ╚═╝ ╚═════╝")
        print("")
        print("--------- K-Effective Eigenvalue Problem ---------")
        print("Outter Solver: Power Iteration")
        print("Inner Sovler:", solver)
        print("Material: ", qmc_data.material_code)
        print("Random Number Generator: ", qmc_data.generator)
        print("Number of Particles per Iteration: ", qmc_data.N)
        print("Number of Spatial Cells: ", qmc_data.Nx)
        print("Initial K: ", qmc_data.keff)

    # iterate over k effective
    while (itt <= max_outter_itt) and (dk >= outter_tol):
        # iterate over scattering source
        phi_new = InnerIteration(qmc_data, solver=solver,
                                 maxit=max_inner_itt, tol=inner_tol,
                                 report_progress=report_progress)
        # phi_hist.append(phi_new)
        k_old = k
        k = UpdateK(phi_old, phi_new, qmc_data)
        k_hist.append(k)
        qmc_data.keff = k
        # res_hist.append(np.linalg.norm(phi_new-phi_old))
        qmc_data.tallies.phi_f = phi_new.copy()
        phi_old = phi_new.copy()  # /norm(phi_new)
        if (qmc_data.source_tilt):
            qmc_data.tallies.dphi_f = qmc_data.tallies.dphi_s
        dk = abs(k - k_old)
        itt += 1
        if (rank == 0) and (report_progress):
            print("**********************")
            print("Iteration:", itt)
            print("k: ", k)
            print("dk: ", dk)
    if (rank == 0):
        if (itt >= max_outter_itt):
            print(
                "Power Iteration convergence to tolerance not achieved: Maximum number of iterations.")
        elif (dk <= outter_tol):
            print("-------------------------------")
            print("Successful Power Iteration convergence.")

    return phi_new, k_hist, itt  # , res_hist


# =============================================================================
# Inner Source Iteration for Power Iteration
# =============================================================================
# TODO: make exitCode an actual output from Picard
def InnerIteration(
        qmc_data,
        solver="LGMRES",
        tol=1e-5,
        maxit=50,
        save_data=False,
        report_progress=True):
    """
    Parameters
    ----------
    qmc_data : TYPE
        DESCRIPTION.
    tol : TYPE, optional
        DESCRIPTION. The default is 1e-5.
    maxit : TYPE, optional
        DESCRIPTION. The default is 50.

    Returns
    -------
    phi : TYPE
        DESCRIPTION.

    """
    comm = MPI.COMM_WORLD
    rank = comm.Get_rank()
    nproc = comm.Get_size()
    Nx = qmc_data.Nx
    G = qmc_data.G
    Nv = Nx * G
    Nt = qmc_data.Nt
    start = time.time()
    matvec_data = MatVec_data(qmc_data)
    if (qmc_data.source_tilt):
        phi0 = np.append(qmc_data.tallies.phi_avg, qmc_data.tallies.dphi_s)
    else:
        phi0 = qmc_data.tallies.phi_avg
    phi0 = np.reshape(phi0, (Nt, 1))

    if (rank == 0) and (report_progress):
        print("     Inner Iteration: ")
    if (solver == "Picard"):
        phi = Picard(qmc_data, tol=tol, maxit=maxit, save_data=False,
                     report_progress=report_progress)
        exitCode = 0
    else:
        A = LinearOperator((Nt, Nt),
                           matvec=MatVec,
                           rmatvec=MatVec,
                           matmat=MatVec,
                           rmatmat=MatVec,
                           dtype=float)
        b = matvec_data[0]
        if (solver == "LGMRES"):
            counter = gmres_counter(disp=report_progress)
            gmres_out = lgmres(
                A,
                b,
                x0=phi0,
                tol=tol,
                maxiter=maxit,
                callback=counter)
        elif (solver == "GMRES"):
            counter = gmres_counter(disp=report_progress)
            gmres_out = gmres(
                A,
                b,
                x0=phi0,
                tol=tol,
                maxiter=maxit,
                callback=counter)
        elif (solver == "BICGSTAB"):
            counter = gmres_counter(disp=report_progress)
            gmres_out = bicgstab(
                A,
                b,
                x0=phi0,
                tol=tol,
                maxiter=maxit,
                callback=counter)
        else:
            print(" Not a valid solver ")
            Exception
        phi = gmres_out[0]
        exitCode = gmres_out[1]
    stop = time.time()
    run_time = stop - start
    if (qmc_data.source_tilt):
        phi = phi[:Nv]
    phi = np.reshape(phi, (Nx, G))

    if (rank == 0):
        if (save_data):
            sim_data = SimData(phi, run_time, tol, nproc)
            SaveData(qmc_data, sim_data)
        if (exitCode > 0) and (report_progress):
            print(
                "     Convergence to tolerance not achieved: Maximum number of iterations.")
        elif (exitCode < 0) and (report_progress):
            print("     Illegal input or breakdown.")
        elif (exitCode == 0) and (report_progress):
            print("     Successful convergence.")

    return phi


def UpdateK(phi_f, phi_s, qmc_data):
    keff = qmc_data.keff
    material = qmc_data.material
    keff *= (np.sum(material.nu * material.sigf * phi_s)
             / np.sum(material.nu * material.sigf * phi_f))
    return keff

# =============================================================================
# Davidson's Algorithm
# =============================================================================
# TODO: Correct normalization of scalar flux in Davidson's output
# TODO: Enable Source Tilting with Davidson's


def Davidson(qmc_data, k0=1.0, l=1, m=None, numSweeps=8, tol=1e-6, maxit=30,
             report_progress=True):
    """
    Parameters
    ----------
    qmc_data : qmc_data structure
    k0 : Float, optional
        DESCRIPTION. The default is 1.0.
    l : Int, optional
        DESCRIPTION. Number of eigenvalues and vectors to solver for The default is 1.
    m : Int, optional
        DESCRIPTION. Restart parameter. The default is 5.
    numSweeps : Int, optional
        DESCRIPTION. The default is 5.
    tol : Float, optional
        DESCRIPTION. The default is 1e-6.
    maxit : Int, optional
        DESCRIPTION. The default is 30.

    Returns
    -------
    phi : TYPE
        DESCRIPTION.
    keff :  TYPE
        DESCRIPTION.
    itt : TYPE
        DESCRIPTION.

    """
    comm = MPI.COMM_WORLD
    rank = comm.Get_rank()

    # Davidson Parameters
    Nt = qmc_data.Nt
    if (qmc_data.source_tilt):
        phi0 = np.append(qmc_data.tallies.phi_avg, qmc_data.tallies.dphi_s)
    else:
        phi0 = qmc_data.tallies.phi_avg
    phi0 = np.reshape(phi0, (Nt))

    # u       = qmc_data.tallies.phi_f.reshape(Nt)
    # orthonormalize initial guess
    V0 = np.array(phi0 / np.linalg.norm(phi0).T)
    V = np.zeros((Nt, maxit))
    axv = np.zeros((Nt, maxit))
    bxv = np.zeros((Nt, maxit))
    Vsize = 1
    V[:, 0] = V0
    k_old = 0.0
    dk = 1.0
    itt = 1
    if (rank == 0):
        print("")
        print("    ██╗ ██████╗ ███╗   ███╗ ██████╗")
        print("    ║ ║██╔═══██╗████╗ ████║██╔════╝")
        print("    ██║██║   ██║██╔████╔██║██║     ")
        print("    ██║██║▄▄ ██║██║╚██╔╝██║██║     ")
        print("    ██║╚██████╔╝██║ ╚═╝ ██║╚██████╗")
        print("    ╚═╝ ╚══▀▀═╝ ╚═╝     ╚═╝ ╚═════╝")
        print("")
        print("--------- K-Effective Eigenvalue Problem ---------")
        print("Outter Solver: Davidson's Method")
        print("Material: ", qmc_data.material_code)
        print("Random Number Generator: ", qmc_data.generator)
        print("Number of Particles per Iteration: ", qmc_data.N)
        print("Number of Spatial Cells: ", qmc_data.Nx)
        print("Initial K: ", qmc_data.keff)

    if (m is None):
        m = maxit + 1  # unless specified there is no restart parameter
    V[:, :Vsize] = PreConditioner(V[:, :Vsize], qmc_data, numSweeps)

    # Davidson Routine
    while (itt <= maxit) and (dk >= tol):
        # print(V)
        if (report_progress):
            print("**********************")
            print(" Davidson Iteration: ", itt)
        axv[:, Vsize - 1] = AxV(V[:, :Vsize], qmc_data)[:, 0]
        bxv[:, Vsize - 1] = BxV(V[:, :Vsize], qmc_data)[:, 0]
        # Scattering linear operator
        AV = np.dot(V[:, :Vsize].T, axv[:, :Vsize])
        BV = np.dot(V[:, :Vsize].T, bxv[:, :Vsize])  # Fission linear operator
        [Lambda, w] = sp.eig(AV, b=BV)  # solve for eigenvalues and vectors
        idx = Lambda.argsort()  # get indices of eigenvalues from smallest to largest
        Lambda = Lambda[idx]       # sort eigenvalues from smalles to largest
        # there can't be any imaginary eigenvalues
        assert (Lambda.imag.all() == 0.0)
        # take the real component of the l largest eigenvalues
        Lambda = Lambda[:l].real
        k = 1 / Lambda
        dk = abs(k - k_old)
        if (report_progress):
            print("K Effective: ", k)
            print("dk: ", dk)
        k_old = k
        w = w[:, idx]          # sort corresponding eigenvector
        w = w[:, :l].real      # take the l largest eigenvectors
        u = np.dot(V[:, :Vsize], w)       # Ritz vectors
        res = AxV(u, qmc_data) - Lambda * BxV(u, qmc_data)  # residual
        t = PreConditioner(res, qmc_data, numSweeps)
        if (Vsize <= m - l):
            Vsize += 1
            # appends new orthogonalization to V
            V[:, :Vsize] = Gram(V[:, :Vsize - 1], t)
        else:
            Vsize = 2
            V[:, :Vsize] = Gram(u, t)  # "restarts" by appending to a new array
        if (itt == maxit):
            print(" Convergence to tolerance not achieved: Maximum number of iterations.")
            break
        else:
            print(" Successful convergence.")
        itt += 1

    keff = 1 / Lambda
    phi = V[:, 0]
    phi = phi / np.linalg.norm(phi).T
    return phi, keff, itt


# =============================================================================
# Functions for Davidson's Method
# =============================================================================

def AxV(V, qmc_data):
    """
    Linear operator for scattering term (I-L^(-1)S)*phi
    """
    v = V[:, -1]
    Nx = qmc_data.Nx
    G = qmc_data.G
    Nt = qmc_data.Nt
    zed = np.zeros((Nx, G))
    phi_in = np.reshape(v, (Nt, 1))
    axv = (phi_in - SI_Map(zed, phi_in, qmc_data))

    return axv


def BxV(V, qmc_data):
    """
    Linear operator for fission term (L^(-1)F*phi)
    """
    v = V[:, -1]
    Nx = qmc_data.Nx
    G = qmc_data.G
    Nv = int(Nx * G)
    Nt = qmc_data.Nt
    zed = np.zeros(Nt)
    phi_in = np.reshape(v, (Nt, 1))
    if (qmc_data.source_tilt):
        dphi = qmc_data.tallies.dphi_s
        qmc_data.tallies.dphi_s = zed

    bxv = SI_Map(phi_in, zed, qmc_data)

    if (qmc_data.source_tilt):
        qmc_data.tallies.dphi_s = dphi
        v[Nv:] = dphi.reshape(Nv)

    return bxv


def PreConditioner(V, qmc_data, numSweeps=8):
    """
    Linear operator approximation of L^(-1)S

    In this case the preconditioner is a specified number of purely scattering
    transport sweeps.
    """
    v = V[:, -1]
    Nx = qmc_data.Nx
    G = qmc_data.G
    Nt = qmc_data.Nt
    Nv = Nx * G
    zed = np.zeros((Nx, G))
    phi_in = np.reshape(v, (Nt, 1))
    for i in range(numSweeps):
        phi_in = SI_Map(zed, phi_in, qmc_data)

    return phi_in


def Gram(V, u):
    """
    Modified Gram Schmidt

    """
    w1 = u - np.dot(V, np.dot(V.T, u))
    v1 = w1 / np.linalg.norm(w1)
    w2 = v1 - np.dot(V, np.dot(V.T, v1))
    v2 = w2 / np.linalg.norm(w2)
    V = np.append(V, v2, axis=1)
    return V

# =============================================================================
# Misc Functions
# =============================================================================


def SimData(phi, time, tol, nproc):
    data = {
        "phi": phi,
        "run_time": time,
        "tolerance": tol,
        "nproc": nproc
    }
    return data
