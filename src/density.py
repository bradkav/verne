"""
density.py - Code for calculating the density of elements in the Earth and atmosphere

Last updated: 17/10/2024
Contact: Bradley Kavanagh, bradkav@gmail.com

"""

import numpy as np
import os
import sys

#--------------------
#Densities are in number/cm^3
#Distances in m
#--------------------

isotopes = None
dens_profiles = None
dens_interp = None
Avals = None
Zvals = None
Niso = None
Niso_full = None
r_list = None

isoID = {"O":0, "Si":1, "Mg":2, "Fe":3, "Ca":4, "Na":5, "S":6, "Al":7, "O_A":8, "N_A": 9}

h_A = 80e3  #Height of atmosphere in (m)
R_E = 6371.0e3  #Earth Radius in (m)

def loadIsotopes():
    print("> VERNE: Loading isotope data and density profiles...")
    
    global dens_profiles, isotopes
    global Avals, Zvals
    global Niso, Niso_full
    global r_list
    
    #Depending on the required directory structure, try 
    #two possible locations for the data files
    rootdir = "data/"
    if not (os.path.exists(rootdir + "isotopes.txt")):
        rootdir = "../data/"
    if not (os.path.exists(rootdir + "isotopes.txt")):
        sys.exit("Data files (isotopes, density profiles, etc.) not found in 'data/' or '../data/'...")
        
    #Load in Earth isotopes
    Zvals = np.loadtxt(rootdir+"isotopes.txt", usecols=(2,))
    Avals = np.loadtxt(rootdir+"isotopes.txt", usecols=(1,)) 
    isotopes = np.loadtxt(rootdir+"isotopes.txt", usecols=(0,))
    Niso = len(isotopes)    #Number of #arth isotopes
    Niso_full = Niso + 2    #Plus isotopes in atmosphere
    
    #Load density profiles for Earth isotopes
    r_list = np.loadtxt(rootdir+"dens_profiles/n_1.dat", usecols=(0,))
    r_list0 = 0.0*r_list
    dens_profiles = [np.loadtxt(rootdir+"dens_profiles/n_"+str(int(iso))+".dat", usecols=(1,)) for iso in isotopes]
    
    #Grid of heights in the atmosphere
    h_list = np.linspace(0, h_A, 100)+1e-5
    
    #Make a slight correction of the profile - truncate at 6371 km
    #This doesn't affect things on large scales, but we have to make
    #sure we've truncated the Earth at R_E exactly if we have very shallow 
    #detectors
    r_list[-1] = R_E
    
    #Append atmospheric points to the list of radii and densities for the Earth
    r_list = np.append(r_list, R_E + h_list)
    for i, dens in enumerate(dens_profiles):
        dens[-1] = dens[-2]
        dens_profiles[i] = np.append(dens, 0.0*h_list)
    
    #Add the atmospheric elements:
    Avals = np.append(Avals,[16, 14])
    Zvals = np.append(Zvals,[8,  7])

    #Load atmospheric parameters and 
    #calculate the atmosphere density profiles...
    Hvals, Tvals, beta = np.loadtxt(rootdir+"ISA.txt", unpack=True)
    Hvals *= 1e3
    beta *= 1e-3 #Get everything in m
    
    #Fraction of Oxygen and Nitrogen
    frac = [0.21, 0.78]
    dens = lambda x: atmos_density(x, Hvals, Tvals, beta)
    dens_vec = np.vectorize(dens)
    dens_atmos = [np.append(r_list0, 2*f_n*dens_vec(h_list)) for f_n in frac]
    dens_profiles.extend(dens_atmos)
    
def dens_interp(index, r):
    return np.interp(r, r_list, dens_profiles[index], left=0, right=0)
    
#International standard atmosphere, ISO 2533:1975
#https://www.iso.org/standard/7472.html
def atmos_density(h, Hvals, Tvals, beta):
    H = R_E*h/(R_E + h)
    if (H > 80000.0):
        return 0
    R = 287.05
    p0 = 1.01e5
    g = 9.807
    
    #determine the layer
    ib = np.digitize(H, Hvals, right=False)-1
    
    if (ib == 0):
        return 0
    if (ib == 1):
        p = p0*(1.0 + beta[1]*H/Tvals[1])**(-g/(beta[1]*R))
    else:
        p = p0*(1.0 + beta[1]*(Hvals[2])/Tvals[1])**(-g/(beta[1]*R))
        for i in range(2,ib):
            if (beta[i] < 1e-3):
                p *= np.exp(-g*(Hvals[i+1]-Hvals[i])/(Tvals[i]*R))
            else:
                p *= (1.0 + beta[i]*(Hvals[i+1]-Hvals[i])/Tvals[i])**(-g/(beta[i]*R))
        if (beta[ib] < 1e-3):
            p *= np.exp(-g*(H-Hvals[ib])/(Tvals[ib]*R))
        else:
            p *= (1.0 + beta[ib]*(H-Hvals[ib])/Tvals[ib])**(-g/(beta[ib]*R))
    n = 1e-6*6.022e23*p/(8314.32e-3*Tvals[ib]) #Air particles per cubic cm
    
    return n
    
loadIsotopes()