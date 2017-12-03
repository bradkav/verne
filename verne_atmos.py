import numpy as np
import WIMpy.DMUtils as DMU
from numpy import pi
from scipy.integrate import quad, dblquad, tplquad
from scipy.interpolate import interp1d, interp2d, griddata, RectBivariateSpline,InterpolatedUnivariateSpline
import os
import sys
import scipy.special
from scipy.special import erf
from skmonaco import mcquad
from scipy.integrate import odeint

from tqdm import tqdm

import mcint

ROOT2 = np.sqrt(2.0)


# Some random functions...

def q(A, E_R):
    mN = 1e3*A*0.9315 #in MeV
    
    return np.sqrt(2*mN*E_R*1e-3) #in MeV

def mu(mX, mA):
    return mX*mA*1.0/(mX+mA)

def ERmax(mX, mA, v):
    return 2*(mu(mX, mA)*v)**2/mA
    
#--------------------

isotopes = None
dens_profiles = None
dens_interp = None
Avals = None
vel_funcs = None
vel_funcs_inv = None
Niso = None
Deff_interp = None
phi_interp = None
corr_interp = None
xmax = np.zeros(8)

Hvals = None
Tvals = None
beta = None

isoID = {"O":0, "N":1}

R_E = 637.1e4

#Densities are in number/cm^3
#Distances in m



def loadIsotopes():
    print " Loading isotope data and density profiles (atmosphere)..."
    
    global dens_profiles
    global isotopes
    global Avals
    global dens_interp
    global Niso
    global corr_interp
    
    global Hvals
    global Tvals
    global beta
    
    Avals = [16, 14]
    
    rootdir = "data/"
    
    Hvals, Tvals, beta = np.loadtxt(rootdir+"ISA/temp.txt", unpack=True)
    Hvals *= 1e3
    beta *= 1e-3 #Get everything in m
    isotopes = np.asarray([0, 1])
    Niso = len(isotopes) - 1
    
    h_list = np.linspace(0, 100000, 1000)
    frac = [0.21, 0.78]
    dens_vec = np.vectorize(density)
    dens_profiles = [2*f_n*dens_vec(h_list) for f_n in frac]
    dens_interp = [interp1d(h_list, dens, bounds_error=False, fill_value = 0) for dens in dens_profiles]
    
    corr_interp = [calcFFcorrection(ID) for ID in range(Niso)]
    
    
def density(h):
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
    
#NB: pi - theta is approximately the angle off the normal for the DM particle hitting the atmosphere...
    
def calcFFcorrection(ID):

    A0 = Avals[ID]
    v_vals = np.linspace(0, 800, 200)

    #print v_vals
    corr_fact = v_vals*0.0
    for i, v in enumerate(v_vals):
            corr_fact[i] = quad(lambda x: 2.0*x*DMU.calcSIFormFactor(x*(1e6/(3e5*3e5))*ERmax(1e10, 0.9315*A0, v), A0), 0, 1)[0]
    corr_fact[0] = 1.0
    return interp1d(v_vals, corr_fact, kind='linear')