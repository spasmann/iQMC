#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon Jun 28 14:08:05 2021

@author: sampasmann
"""

import numpy as np
import math
import random



def cycles(Geo,N,Nr,R, Sig_t,Sig_f,Sig_s,Sig_c,X,nu,weights = [],
                       inactive_cycles = 5, active_cycles = 20):
    #total # of cycles
    cycles = inactive_cycles + active_cycles
    #how many groups
    G = np.shape(Sig_t)[1]
    #how many shells
    num_regions = len(R)
    #how many cells
    N = int(N)
    #scalar flux
    scalar_flux = np.zeros((cycles,G,Nr))
    #shannon entropy
    entropy = np.zeros((cycles,Nr))
    
    #macroscopic cross sections (cm-1)
    Sig_a = Sig_c + Sig_f
    Sig_st = np.zeros((num_regions,G))
    
    #total scattering cross section per energy group
    for j in range(num_regions):
        for i in range(G):
            Sig_st[j,i] = np.sum(Sig_s[j,i,:])
            
    #stratified sampling and scalar flux discretization
    dr = R[-1]/Nr
    lowR = np.linspace(0,(R[-1]-dr),Nr)
    highR = np.linspace(dr,(R[-1]),Nr)
    fission_sites = []

    #initial fission sites
    for cell in range(Nr):
        fission_sites = np.append(fission_sites,
                                  np.random.uniform(lowR[cell],highR[cell],N))
    print('Expectation error: ',np.abs(R[-1]/2 - np.mean(fission_sites)))
    
    #total number of simulated particles
    L = int(np.ceil(N/len(fission_sites)))*len(fission_sites)
    
    if (len(weights) == 0):
        weights = np.ones(L)
    else:
        weight_coeff = weights
        weights = []
        assert (len(weight_coeff) == Nr)
        for i in range(Nr):
            weights = np.append(weights,np.ones(N)*weight_coeff[i])
    assert (len(weights) == len(fission_sites))
    
    #radii of fission sites
    radii = fission_sites.copy()
    old_gen = np.sum(weights)
    k = np.zeros(cycles)
    
    #volume of cell for scalar flux calc
    Vcell = lambda p: dr
    if (Geo == 2):
        Vcell = lambda p: math.pi*(highR[p]**2 - lowR[p]**2)
    elif (Geo == 3):
        Vcell = lambda p: (4.0/3.0)*math.pi*(highR[p]**3 - lowR[p]**3)
    
    for cycle in range(cycles):
        print("*********** Cycle ",cycle, " ***********")
        fission_sites = []
        fission_site_weights = []
        #shannon entropy - source sites per element
        if (cycle == 0):
            entropy[cycle,:] = N/(Nr*N)
        else:
            for i in range(len(radii)):
                #source sites in ith element
                cell = (math.fabs(radii[i])>=lowR)*(math.fabs(radii[i])<highR)
                entropy[cycle,cell] += 1
            #divided by total number of source sites
            entropy[cycle,:] /= len(radii)
        
        #assert(weights.size == radii.size)
        #print(cycle, L)
        for neut in range(L):
            #define starting region 
            region = 0
            #define starting group
            group_prob = random.random()
            Group = 0
            X_sum = X[region,Group]
            if (group_prob > X_sum):
                while (group_prob > X_sum):
                    Group += 1
                    X_sum += X[region,Group]
            #grab neutron from stack
            r = radii[neut]
            #which cell is the neutron born in
            cell_i = (math.fabs(r)>=lowR)*(math.fabs(r)<highR)
            #neutron weight
            weight = float(weights[neut])
            #initialize coordinates
            x = y = s = 0.0
            z = r
            #define angles based on geometrey
            #SLAB
            mu = np.random.uniform(-1,1)
            if (Geo > 1):
                #CYLINDER & SPHERE
                phi = np.random.uniform(0,(2*np.pi))
                muSin = math.sqrt(1-mu**2)
                
            #neutron life cycle loop
            alive = 1
            while (alive):
                #get distance to collision
                s = -math.log(1.0-random.random())/Sig_t[region,Group]
                distance_to_edge = s 
                collide = 1          
                #check distance to boundaries
                #Slab
                if (Geo == 1):
                    for w in range(region,region+1):
                        distance_to_edge = (R[w]-r)/mu
                        if (distance_to_edge > 0) and (distance_to_edge < s):
                            s = distance_to_edge + 1e-10
                            collide = 0
                    #center reflection
                    if (region == 0):
                        #determine region
                        distance_to_center = (0-r)/mu
                        if (0 < distance_to_center < s):
                            collide = 0
                            s = -distance_to_center + 1e-10
                            mu *= -1
                else:
                    #Cylinder
                    if (Geo == 2):
                        a = muSin*math.cos(phi)**2 + muSin*math.sin(phi)**2
                        h = (z*muSin*math.cos(phi) + y*muSin*math.sin(phi))
                    #Sphere
                    else:
                        h = (x*muSin*math.cos(phi) + y*muSin*math.sin(phi) + z*mu)
                        a = 1
                    for w in range(region,region+1):
                        c = (x**2 + y**2 + z**2 - R[w]**2)
                        if (a != 0) and (h**2 - a*c > 0):
                            d1 = (-h + np.sqrt(h**2 - a*c))/a
                            d2 = (-h - np.sqrt(h**2 - a*c))/a
                            if (c < 0):
                                if (d1 > 0):
                                    distance_to_edge = d1
                                else:
                                    distance_to_edge = d2
                            elif (d1 > 0):
                                if (d1 > d2):
                                    distance_to_edge = d2
                                else:
                                    distance_to_edge = d1
                            if (distance_to_edge < s):
                                s = distance_to_edge + 1e-10
                                collide = 0
                #update position
                #CYLINDER
                if (Geo == 2):
                    z += s*muSin*math.cos(phi)
                    y += s*muSin*math.sin(phi)
                #SLAB
                else:
                    z += s*mu
                    #SPHERE
                    if (Geo == 3):
                        x += s*muSin*math.cos(phi)
                        y += s*muSin*math.sin(phi)
                    
                #move particle
                r = math.sqrt(x**2 + y**2 + z**2)
                
                #determine region
                if (r < R[0]):
                    region = 0
                else:
                    for w in range(num_regions-1):
                        if (R[w] < r < R[w+1]):
                            region = w+1
                            
                #still in the shell?
                if (math.fabs(r) >= R[num_regions-1]): 
                    alive = 0
                    
                elif (collide):
                    #check which cell the neutron is in
                    cell_f = (math.fabs(r)>=lowR)*(math.fabs(r)<highR)
                    #calculate scalar flux                        
                    scalar_flux[cycle,Group,cell_f] += (weight/
                                            (Sig_st[region,Group]*Vcell(cell_f)))
                    #decide if collision is abs or scat
                    scatter_prob = random.random()
                    if (scatter_prob < Sig_st[region,Group]/Sig_t[region,Group]):
                        #pick new angles based on geometrey
                        #SLAB
                        mu = np.random.uniform(-1,1)
                        if (Geo > 1):
                            #CYLINDER & SPHERE
                            phi = np.random.uniform(0,(2*np.pi))
                            muSin = math.sqrt(1-mu**2)
                        #pick group
                        scatter_prob = random.random()
                        scatter_sum = Sig_s[region,Group,0]/Sig_st[region,Group]
                        if (scatter_prob < scatter_sum):
                            Group = 0
                        else:
                            i = 1
                            while (scatter_prob > scatter_sum):
                                scatter_sum += (Sig_s[region,Group,i]/
                                                         Sig_st[region,Group])
                                scatter_group = i
                                i += 1                                
                            Group = scatter_group
                            #weight = weights[neut]*nu[region,Group]
                    else:
                        fiss_prob = random.random()
                        alive = 0
                        if (fiss_prob <= Sig_f[region,Group]/Sig_a[region,Group]):
                            #fission
                            fission_sites.append(r)
                            weight *= nu[region,Group]
                            fission_site_weights.append(weight)
                                
                            
        #sample neutrons for next generation from fission sites
        num_per_site = int(np.ceil(N*Nr/len(fission_sites)))
        #delete the initial sites and weights
        radii = []
        weights = []
        L = num_per_site*len(fission_sites)
        #replace radii and weights with new values
        for site in range(len(fission_sites)):
            for addl in range(num_per_site):
                radii.append(fission_sites[site])
                weights.append(fission_site_weights[site]/num_per_site)
        
        radii =  np.array(radii)

            
        weights =  np.array(weights)
        new_gen = np.sum(weights)
        #k is ratio of neutrons per generation
        k[cycle] = new_gen/old_gen
        old_gen = new_gen
    
    
    keff = np.mean(k[inactive_cycles:])
    scalar_flux = np.mean(scalar_flux[inactive_cycles:,:,:],axis=0)
    scalar_flux /= scalar_flux[0,0]
    
    #shannon entropy - 
    E = np.zeros(cycles)
    for cycle in range(cycles):
        for s in range(Nr):
            E[cycle] -= entropy[cycle,s]*np.log2(entropy[cycle,s])

    return keff, scalar_flux, E
