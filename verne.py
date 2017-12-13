import numpy as np
from scipy.integrate import quad, simps
from scipy.interpolate import interp1d, interp2d
from scipy.integrate import odeint
import scipy.special
import os
import sys

import MaxwellBoltzmann as MB
    
#--------------------
#Theta = 0 is directly from BELOW, angles in radians
    
#Densities are in number/cm^3
#Distances in m
#--------------------

isotopes = None
dens_profiles = None
dens_interp = None
Avals = None
Niso = None
Niso_full = None

phi_interp = None
corr_interp = None
corr_Pb = None
corr_Cu = None

isoID = {"O":0, "Si":1, "Mg":2, "Fe":3, "Ca":4, "Na":5, "S":6, "Al":7, "O_A":8, "N_A": 9}

h_A = 80e3  #Height of atmosphere in (m)
R_E = 6371.0e3  #Earth Radius in (m)

def loadIsotopes():
    print "    Loading isotope data and density profiles..."
    
    global dens_profiles
    global isotopes
    global Avals
    global dens_interp
    global Niso
    global Niso_full
    
    rootdir = "data/"
    
    #Load in Earth isotopes
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
    #sure we've truncated the Earth at R_E if we have very shallow 
    #detectors
    r_list[-1] = R_E
    
    #Append atmospheric points to the list of radii and densities for the Earth
    r_list = np.append(r_list, R_E + h_list)
    for i, dens in enumerate(dens_profiles):
        dens[-1] = dens[-2]
        dens_profiles[i] = np.append(dens, 0.0*h_list)
    
    #Add the atmospheric elements:
    Avals = np.append(Avals,[16, 14])

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
    
    #Generate interpolation functions for the density profiles
    dens_interp = [interp1d(r_list, dens, bounds_error=False, fill_value = 0) for dens in dens_profiles]
    
#Generate interpolation functions for the Form Factor corrections (C_i(v))
def loadFFcorrections(m_x):
    global corr_interp
    global corr_Pb
    global corr_Cu
    
    #Check that the isotope list has been loaded
    if (Avals is None):
        loadIsotopes()
    
    print "    Calculating Form Factor corrections for m_x = ", m_x, " GeV..."
    corr_interp = [calcFFcorrection(m_x, Avals[ID]) for ID in range(Niso_full)]
    corr_Pb = calcFFcorrection(m_x, 207) #Also need Lead + Copper, for the shielding
    corr_Cu = calcFFcorrection(m_x, 63.5)
    
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
    
    
#Path length [m], as measured from the top of the atmosphere to the detector 
#(at 'depth' m underground)
def pathLength(depth, theta):
    r_det = R_E - depth
    return +np.cos(theta)*r_det + np.sqrt((-np.cos(theta)*r_det)**2 - (r_det**2 - (R_E+h_A)**2))
   
#Path length [m], as measured from the Earth's surface to the detector 
#(at 'depth' m underground)
def pathLength_Earth(depth, theta):
    r_det = R_E - depth
    return +np.cos(theta)*r_det + np.sqrt((-np.cos(theta)*r_det)**2 - (r_det**2 - (R_E)**2))
    

#Calculate the Form Factor correction for a nucleus of mass-number A0
#See BJK:give equation number


#Maximum recoil energy (in keV)
def ERmax(mX, mA, v):
    mu = mX*mA*1.0/(mX+mA)
    return (1e6/(3e5*3e5))*2*(mu*v)**2/mA

#Calculate the spin-independent form factor 
#for nucleon number A0 and recoil energy E
def calcSIFormFactor(E, A0):
        #Helm
        if (E < 1e-5):
            return 1.0

        #Define conversion factor from amu-->keV
        amu = 931.5*1e3

        #Convert recoil energy to momentum transfer q in keV
        q1 = np.sqrt(2*A0*amu*E)

        #Convert q into fm^-1
        q2 = q1*(1e-12/1.97e-7)
    
        #Calculate nuclear parameters
        s = 0.9
        a = 0.52
        c = 1.23*(A0**(1.0/3.0)) - 0.60
        R1 = np.sqrt(c*c + 7*np.pi*np.pi*a*a/3.0 - 5*s*s)
    
        x = q2*R1
        J1 = np.sin(x)/x**2 - np.cos(x)/x
        F = 3*J1/x
        return (F**2)*(np.exp(-(q2*s)**2))

def calcFFcorrection(m_x, A0):
    v_vals = np.linspace(0, 1000, 200)
    corr_fact = v_vals*0.0
    for i, v in enumerate(v_vals):
        corr_fact[i] = quad(lambda x: 2.0*x*calcSIFormFactor(x*ERmax(m_x, 0.9315*A0, v), A0), 0, 1)[0]
    corr_fact[0] = 1.0
    return interp1d(v_vals, corr_fact, kind='linear', bounds_error=False, fill_value=0.0)

#Calculate the DM-nucleus 'effective' cross section
#which takes into account the average energy loss
def effectiveXS(sigma_p, m_X, A):
    m_p = 0.9315 #Proton mass
    m_A = 0.9315*A
    mu_A = m_A*m_X/(m_A + m_X)
    mu_p = m_p*m_X/(m_p + m_X)
    
    return sigma_p*(1.0/(m_X*m_A))*(A**2)*(mu_A**4/mu_p**2)
    

def CalcF(vf, gamma, depth,sigma_p, m_x, target, vmax_interp):
    
    #Define a grid of values for theta which we sample over
    #theta = pi/2 is often problematic, so we sample more densely there
    tlist = np.linspace(0, np.pi, 101)
    tlist = np.append(tlist, (np.pi/2)*(1 + np.logspace(-3, -0.01, 50)))
    tlist = np.append(tlist, (np.pi/2)*(1 - np.logspace(-3, -0.01, 50)))
    tlist = np.sort(tlist)
    
    fint = tlist*0.0
    for i in range(len(tlist)):
        #If maximum vf you can get this value of theta is greater
        #than the speed we're interested in, set to zero
        if (vmax_interp(tlist[i]) < vf): 
            fint[i] = 0.0
        else:
            fint[i] = f_integrand_full(vf, tlist[i], gamma, depth, sigma_p, m_x, target)

    #Integrate with Simpson's rule
    return simps(fint, tlist)
    
    
    
def f_integrand_full(vf, theta, gamma, depth, sigma_p, m_x, target):

    #Calculate the initial velocity corresponding to this final velocity vf
    dv = 1.5
    vi1 = calcVinitial_full(vf+dv/2.0, theta,  depth, sigma_p, m_x, target)
    vi2 = calcVinitial_full(vf-dv/2.0, theta,  depth, sigma_p, m_x, target)

    #Calculate the average and the numerical derivative
    vi = (vi1 + vi2)/2.0
    dvi_by_dvf = np.abs(vi1 - vi2)*1.0/dv

    return (dvi_by_dvf)*np.sin(theta)*(vi**2)*MB.calcf_integ(vi, theta, gamma)
 
#Calculate the distance of a point from the centre of the Earth
#The point is defined by:
#   - theta, the angle of the trajectory
#   - depth,the detector depth
#   - D, the distance along the trajectory, starting at the top of the atmosphere

def radius(D, theta, depth):
    r_det = R_E - depth 
    return np.sqrt((R_E+h_A)**2 + D**2 + 2*D*(r_det*np.cos(theta) - pathLength(depth, theta)))

#Derivative of DM speed along path length D
#To be used by the ODE integrator
def dv_by_dD(v, D, params):

    theta, depth, sigma_p, m_x, target = params
    res = 0.0
    if (target == "atmos"):
        isovals = [8,9]
    elif (target == "earth"):
        isovals = range(Niso)
    else:
        isovals = range(Niso_full)
    
    r = radius(D, theta, depth)

    #Loop over the relevant isotops
    for i in isovals:
        res += dens_interp[i](r)*effectiveXS(sigma_p, m_x, Avals[i])*corr_interp[i](v)
    return -1e2*v*res #(km/s)/m


def dv_by_dD_Pb(v, D, params):        
    #Pb density
    n_Pb = 3.3e22
    A_Pb = 207
        
    sigma_p, m_x = params
    res = n_Pb*effectiveXS(sigma_p, m_x, A_Pb)*corr_Pb(v)
    return -1e2*v*res #(km/s)/m
    
def dv_by_dD_Cu(v, D, params):        
    #Pb density
    #n_Pb = 3.3e22
    n_Cu = 8.5e22
    A_Cu = 63.5
        
    sigma_p, m_x = params
    res = n_Cu*effectiveXS(sigma_p, m_x, A_Cu)*corr_Cu(v)
    return -1e2*v*res #(km/s)/m

def calcVfinal(v0, theta,  depth, sigma_p, m_x, target="full"):
    params = [theta, depth, sigma_p, m_x, target]

    #Propagate across the atmosphere
    if (target == "atmos"):
        d1 = 0
        d2 = pathLength(depth, theta) - pathLength_Earth(depth, theta)

    #Propagate from the surface of the Earth to the detector
    if (target == "earth"):
        d1 = pathLength(depth, theta) - pathLength_Earth(depth, theta)
        d2 = pathLength(depth, theta)
    
    psoln = odeint(dv_by_dD, v0, [d1,d2] , args=(params,), mxstep=1000, rtol=1e-6)
    vf = psoln[1]
    return vf
    
def calcVfinal_full(v0, theta,  depth, sigma_p, m_x, target="full"):
    vf = 1.0*v0
    if (target in ["atmos", "full", "no_shield", "SUF", "MPI"]):
        vf = calcVfinal(vf, theta,  depth, sigma_p, m_x, target="atmos")
    if (target in ["earth", "full", "no_shield", "SUF", "MPI"]):
        vf = calcVfinal(vf, theta,  depth, sigma_p, m_x, target="earth") 
    if (target == "MPI"):
        vf = calcVfinal_shield_MPI(vf, sigma_p, m_x)
    if (target == "SUF"):
        vf = calcVfinal_shield_SUF(vf, sigma_p, m_x)
    return vf
    
def calcVinitial(v0, theta,  depth, sigma_p, m_x, target="earth"):
    params = [theta, depth, sigma_p, m_x, target]

    #Propagate across the atmosphere
    if (target == "atmos"):
        d1 = 0
        d2 = pathLength(depth, theta) - pathLength_Earth(depth, theta)

    #Propagate from the surface of the Earth to the detector
    if (target == "earth"):
        d1 = pathLength(depth, theta) - pathLength_Earth(depth, theta)
        d2 = pathLength(depth, theta)

    
    psoln = odeint(dv_by_dD, v0, [d2,d1], args=(params,), mxstep=1000, rtol=1e-6)
    #print v0, psoln[1]
    return psoln[1]
    #Initial values
    
def calcVinitial_full(v0, theta,  depth, sigma_p, m_x, target="full"):
    vi = 1.0*v0
    if (target == "MPI"):
        vi = calcVinitial_shield_MPI(vi, sigma_p, m_x)
    if (target == "SUF"):
        vi = calcVinitial_shield_SUF(vi, sigma_p, m_x)
    if (target in ["earth", "full", "no_shield", "SUF", "MPI"]):
        vi = calcVinitial(vi, theta,  depth, sigma_p, m_x, target="earth")
    if (target in ["atmos", "full", "no_shield", "SUF", "MPI"]):
        vi = calcVinitial(vi, theta,  depth, sigma_p, m_x, target="atmos")

    return vi
    
def calcVfinal_shield_SUF(v0, sigma_p, m_x):
    params = [sigma_p, m_x]
    #Propagate through 16cm of Lead
    psoln = odeint(dv_by_dD_Pb, v0, [0,16.0e-2] , args=(params,), mxstep=1000, rtol=1e-6)
    return psoln[1]
    
def calcVinitial_shield_SUF(v0,  sigma_p, m_x):
    params = [sigma_p, m_x]
    #Propagate through 16cm of Lead (backwards)
    psoln = odeint(dv_by_dD_Pb, v0, [16.0e-2,0] , args=(params,), mxstep=1000, rtol=1e-6)
    return psoln[1]

def calcVfinal_shield_MPI(v0, sigma_p, m_x):
    params = [sigma_p, m_x]
    #Propagate through 1mm Copper
    psoln = odeint(dv_by_dD_Cu, v0, [0,1e-3] , args=(params,), mxstep=1000, rtol=1e-6)
    return psoln[1]
    
def calcVinitial_shield_MPI(v0, sigma_p, m_x):
    params = [sigma_p, m_x]
    #Propagate through 1mm Copper
    psoln = odeint(dv_by_dD_Cu, v0, [1e-3,0] , args=(params,), mxstep=1000, rtol=1e-6)
    return psoln[1]

    