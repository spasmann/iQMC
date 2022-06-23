#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed May 11 12:08:53 2022

@author: sampasmann
"""

import sys
sys.path.append("../")
import os
from src.init_files.mg_init import MultiGroupInit
import numpy as np
import matplotlib.pyplot as plt

Nx = 1
data12 = MultiGroupInit(numGroups=12, Nx=Nx)
data70 = MultiGroupInit(numGroups=70, Nx=Nx)
data618 = MultiGroupInit(numGroups=618, Nx=Nx)

script_dir = os.path.dirname(__file__)
rel_path = "../src/materials/HDPE/"
abs_file_path = os.path.join(script_dir, rel_path)
centers12 = np.genfromtxt(abs_file_path+"group_centers_12G_HDPE.csv", delimiter=",")
centers70 = np.genfromtxt(abs_file_path+"group_centers_70G_HDPE.csv", delimiter=",")
centers618 = np.genfromtxt(abs_file_path+"group_centers_618G_HDPE.csv", delimiter=",")

edges12 = np.genfromtxt(abs_file_path+"group_edges_12G_HDPE.csv", delimiter=",")
edges70 = np.genfromtxt(abs_file_path+"group_edges_70G_HDPE.csv", delimiter=",")
edges618 = np.genfromtxt(abs_file_path+"group_edges_618G_HDPE.csv", delimiter=",")

dE12 = abs(edges12[1:] - edges12[:-1])
dE70 = abs(edges70[1:] - edges70[:-1])
dE618 = abs(edges618[1:] - edges618[:-1])

y12 = (data12.true_flux[0,:]/dE12)
y12 /= np.sum(y12*dE12)
y70 = (data70.true_flux[0,:]/dE70)
y70 /= np.sum(y70*dE70)
y618 = (data618.true_flux[0,:]/dE618)
y618 /= np.sum(y618*dE618)

centers12[-1] = np.min(centers618)
centers70[-1] = np.min(centers618)

plt.figure(dpi=300)
size = 3
where = 'mid'
drawstyle='steps-mid'
plt.step(centers12, y12, '-o', where=where,drawstyle=drawstyle, markersize=size,label='G=12')
plt.step(centers70, y70, where=where, drawstyle=drawstyle, markersize=size,label='G=70')
plt.step(centers618, y618,where=where, drawstyle=drawstyle, markersize=size,label='G=618')

plt.legend()
plt.xscale('log')
plt.yscale('log')
plt.xlabel('E (MeV)')
plt.ylabel(r'$\phi(E)$')

