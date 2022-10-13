"""
verne.py - Code for calculating the Earth stopping effect, primarily for heavy Dark Matter.


Last updated: 13/12/2017
Contact: Bradley Kavanagh, bradkav@gmail.com

"""

import numpy as np
from scipy.integrate import quad, simps
from scipy.interpolate import interp1d, interp2d
from scipy.integrate import odeint

import scipy.special
import os
import sys

import MaxwellBoltzmann as MB

from scipy.integrate import solve_ivp

def odeint_new(f, v0, t_span , args, mxstep, rtol):
    #print(v0)
    result = solve_ivp(f, t_span, y0=np.atleast_1d(v0), args=args, rtol=rtol, method="RK23")
    #print(result)
    #print(result.y[0][-1])
    return [result.t[-1], result.y[0][-1]]


#--------------------
#Theta = 0 is directly from BELOW, angles in radians
#Gamma = 0 is mean DM flux directly from BELOW, angles in radians

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

phi_interp = None

corr_interp = None
corr_Pb = None
corr_Cu = None
corr_Fe = None

#Note that Form Factors are only implemented for spin-independent interactions at the moment
NEGLECT_FF = True
if (NEGLECT_FF == True):
    print("> verne.py: Neglecting form factors...")

isoID = {"O":0, "Si":1, "Mg":2, "Fe":3, "Ca":4, "Na":5, "S":6, "Al":7, "O_A":8, "N_A": 9}

h_A = 80e3  #Height of atmosphere in (m)
R_E = 6371.0e3  #Earth Radius in (m)

m_e = 511e-3 #Electron mass in GeV
alpha = 1/137.03 #Fine structure constant
q_screen = alpha*m_e #Electron screening momentum in GeV

#--------------------
#Integration parameters

TOL = 1e-5
NSTEP = 200


#--------------------


def loadIsotopes():
    print("> VERNE: Loading isotope data and density profiles...")
    
    global dens_profiles
    global isotopes
    global Avals, Zvals
    global dens_interp
    global Niso
    global Niso_full
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
    #Generate interpolation functions for the density profiles
    #dens_interp = [interp1d(r_list, dens, bounds_error=False, fill_value = 0) for dens in dens_profiles]
    
def dens_interp(index, r):
    return np.interp(r, r_list, dens_profiles[index], left=0, right=0)
    
    
#Generate interpolation functions for the Form Factor corrections (C_i(v))
def loadFFcorrections(m_x, interaction = "SI"):
    global corr_interp, corr_Pb, corr_Cu, corr_Fe

   
    #Check that the isotope list has been loaded
    if (Avals is None):
        loadIsotopes()
    
    print("> VERNE: Calculating Form Factor corrections for m_x = ", m_x, " GeV, with " + interaction + " interactions...")
    corr_interp = [calcFFcorrection(m_x, Avals[ID], interaction) for ID in range(Niso_full)]
    #Also need Lead + Copper, for the shielding
    corr_Pb = calcFFcorrection(m_x, 207, interaction) 
    corr_Cu = calcFFcorrection(m_x, 63.5, interaction)
    corr_Fe = calcFFcorrection(m_x, 55.8, interaction)

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
#See Eq. (10) of the paper

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

print("Add hDP...")
def FFcorrection_integrand(x, v, m_x, A0, interaction="SI"):
    if (interaction.lower() == "SI".lower() or interaction.lower() == "hDP".lower()):
        return 2.0*x*calcSIFormFactor(x*ERmax(m_x, 0.9315*A0, v), A0)
        
    elif (interaction.lower() == "Millicharge".lower()):
        m_A = 0.9315*A0
        mu_A = m_A*m_x/(m_A + m_x)
        q_max = 2*mu_A*v/3e5
        x_s = (q_screen/q_max)**2
        return (1/np.log(1/x_s))*(1/x)*calcSIFormFactor(x*ERmax(m_x, 0.9315*A0, v), A0)
        
    else:
        return 1.0

def calcFFcorrection(m_x, A0, interaction = "SI"):
    if (NEGLECT_FF):
        return lambda x: 1.0
    
    v_vals = np.linspace(0.001, 1000, 200)
    corr_fact = v_vals*0.0
    
    x_min = 0
    x_max = 1
    for i, v in enumerate(v_vals):
        if (interaction.lower() == "Millicharge".lower()):
            m_A = 0.9315*A0
            mu_A = m_A*m_x/(m_A + m_x)
            q_max = 2*mu_A*v/3e5
            x_s = (q_screen/q_max)**2
            x_min = x_s
            x_max = 1
        
        if (x_min < x_max):
            corr_fact[i] = quad(lambda x: FFcorrection_integrand(x, v, m_x, A0, interaction), x_min, x_max)[0]
        else:
            corr_fact[i] = 0.0
            
    if (interaction.lower() == "SI".lower() or interaction.lower() == "hDP".lower()):
        corr_fact[0] = 1.0
    elif (interaction.lower() == "Millicharge".lower()):
        corr_fact[0] = 0.0
        
    
    return interp1d(v_vals, corr_fact, kind='linear', bounds_error=False, fill_value=0.0)


#Calculate the DM-nucleus 'effective' cross section
#which takes into account the average energy loss
#  sigma_effective = sigma*<E_R>/(m_x v^2)
def effectiveXS(sigma_p, m_X, A, Z, v, interaction="SI"):
    m_p = 0.9315 #Proton mass
    m_A = m_p*A
    mu_A = m_A*m_X/(m_A + m_X)
    mu_p = m_p*m_X/(m_p + m_X)
    
    #C is 'interaction enhancement' factor
    if (interaction.lower() == "SI".lower()):
        C = A**2
        
    elif (interaction.lower() == "hDP".lower()):
        C = Z**2 
        
    elif (interaction.lower() == "SD".lower()):
        #Include SD only for nitrogen - valid for coupling purely to protons OR neutrons only
        #Neglecting a subdominant contribution from Oxygen-17 in the atmosphere
        if (A == 14):
            J_N = 1.0
            S = 0.5
            C = (4./3.)*((J_N + 1)/J_N)*S**2
        else:
            C = 0.0
    elif (interaction.lower() == "Millicharge".lower()):
        v_s = 3e5*alpha*m_e/(2*mu_A) #Screening velocity
        if (v < v_s):
            return 0.0
        C = Z**2*(2*v_s/v)**4*np.log(v/v_s)
    else:
        raise ValueError("Cross sections only defined for interactions: 'SI', 'SD', 'Millicharge'...")
    
    return sigma_p*(1.0/(m_X*m_A))*C*(mu_A**4/mu_p**2)
    
    

#Calculate the final speed distribution at the detector
def CalcF(vf, gamma, depth,sigma_p, m_x, target, vmax_interp, interaction="SI"):
    
    #Define a grid of values for theta which we sample over
    #theta = pi/2 is often problematic, so we sample more densely there
    tlist = np.linspace(0, np.pi, 101)
    tlist = np.append(tlist, (np.pi/2)*(1 + np.logspace(-3, -0.01, 50)))
    tlist = np.append(tlist, (np.pi/2)*(1 - np.logspace(-3, -0.01, 50)))
    tlist = np.sort(tlist)
    
    fint = tlist*0.0
    for i in range(len(tlist)):
        #If maximum vf you can get for this value of theta is greater
        #than the speed we're interested in, set to zero
        if (vmax_interp(tlist[i]) < vf): 
            fint[i] = 0.0
        else:
            fint[i] = f_integrand_full(vf, tlist[i], gamma, depth, sigma_p, m_x, interaction, target)

    #Integrate with Simpson's rule
    return simps(fint, tlist)
    
    
#Integrand for calculating the final speed distribution at the detector
def f_integrand_full(vf, theta, gamma, depth, sigma_p, m_x, interaction, target):
    #Calculate the initial velocity corresponding to this final velocity vf
    dv = 0.5
    vi1 = calcVinitial_full(vf-dv/2.0, theta,  depth, sigma_p, m_x, interaction, target)
    vi2 = calcVinitial_full(vf+dv/2.0, theta,  depth, sigma_p, m_x, interaction, target)
    #vi3 = calcVinitial_full(vf+dv/2.0, theta,  depth, sigma_p, m_x, interaction, target)
    #vi4 = calcVinitial_full(vf+dv, theta,  depth, sigma_p, m_x, interaction, target)    

    #Calculate the average and the numerical derivative
    vi = (vi1 + vi2)/2.0
    dvi_by_dvf = np.abs(vi1 - vi2)*1.0/dv
    #vi = (vi1 + vi2 + vi3 + vi4)/4.0
    #dvi_by_dvf = np.abs(-vi4 + 8*vi3 - 8*vi2 + vi1)*1.0/(6*dv)
    
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

    theta, depth, sigma_p, m_x, interaction, target = params
    res = 0.0
    if (target == "atmos"):
        isovals = [8,9]
    elif (target == "earth"):
        isovals = range(Niso)
    else:
        isovals = range(Niso_full)
    
    r = radius(D, theta, depth)

    #Loop over the relevant isotopes
    for i in isovals:
        #Only include a relevant form factor correction for spin-independent interactions
        #If another interaction is added which needs a form factor, add the correction here!
        if (interaction.lower() in ["SI".lower(), "hDP".lower(), "Millicharge".lower()]):
            FF_correction = corr_interp[i](v)
        else:
            FF_correction = 1.0
        
        res += dens_interp(i, r)*effectiveXS(sigma_p, m_x, Avals[i], Zvals[i], v, interaction)*FF_correction
    return -1e2*v*res #(km/s)/m

def dv_by_dD_concrete(v, D, params):

    sigma_p, m_x, interaction = params
    isovals = range(Niso)

    #Density of isotopes just below the surface
    r = R_E - 1.0

    res = 0.0
    #Loop over the relevant isotopes                              
    for i in isovals:
        #print(i, dens_interp(i, r))
        if (interaction.lower() in ["SI".lower(), "hDP".lower(), "Millicharge".lower()]):
            FF_correction = corr_interp[i](v)
        else:
            FF_correction = 1.0
        
        res += dens_interp(i, r)*effectiveXS(sigma_p, m_x, Avals[i], Zvals[i], v, interaction)*FF_correction
    return -1e2*v*res #(km/s)/m   

#Derivative for the case of Pb shielding
def dv_by_dD_Pb(v, D, params):        
    #Pb density
    n_Pb = 3.3e22
    A_Pb = 207
    Z_Pb = 82
        
    sigma_p, m_x, interaction = params
    
    if (interaction.lower() in ["SI".lower(), "hDP".lower(), "Millicharge".lower()]):
        FF_correction = corr_Pb(v)
    else:
        FF_correction = 1.0

        
    res = n_Pb*effectiveXS(sigma_p, m_x, A_Pb, Z_Pb, v, interaction)*FF_correction
    return -1e2*v*res #(km/s)/m
    
#Derivative for the case of Cu shielding
def dv_by_dD_Cu(v, D, params):        
    #Cu density
    n_Cu = 8.5e22
    A_Cu = 63.5
    Z_Cu = 29
        
    sigma_p, m_x, interaction = params
    
    if (interaction.lower() in ["SI".lower(), "hDP".lower(), "Millicharge".lower()]):
        FF_correction = corr_Cu(v)
    else:
        FF_correction = 1.0
    
    res = n_Cu*effectiveXS(sigma_p, m_x, A_Cu, Z_Cu, v, interaction)*FF_correction
    return -1e2*v*res #(km/s)/m

#Derivative for the case of Fe shielding  
def dv_by_dD_Fe(v, D, params):
    #Fe density                                                                                                                                              
    #n_Fe = 7.874 #g/cm^3
    n_Fe = 8.5e22
    A_Fe = 55.85
    Z_Fe = 26

    sigma_p, m_x, interaction = params
    if (interaction.lower() in ["SI".lower(), "hDP".lower(), "Millicharge".lower()]):
        FF_correction = corr_Fe(v)
    else:
        FF_correction = 1.0
    
    res = n_Fe*effectiveXS(sigma_p, m_x, A_Fe, Z_Fe, v, interaction)*FF_correction
    return -1e2*v*res #(km/s)/m  


#Calculate the final velocity after propagating across 'target'
#Here, target = "atmos" or "earth"
def calcVfinal(vi, theta,  depth, sigma_p, m_x, interaction="SI", target="full"):
    params = [theta, depth, sigma_p, m_x, interaction, target]

    #Propagate across the atmosphere
    if (target == "atmos"):
        d1 = 0
        d2 = pathLength(depth, theta) - pathLength_Earth(depth, theta)

    #Propagate from the surface of the Earth to the detector
    if (target == "earth"):
        d1 = pathLength(depth, theta) - pathLength_Earth(depth, theta)
        d2 = pathLength(depth, theta)
    
    psoln = odeint(dv_by_dD, vi, [d1,d2] , args=(params,), mxstep=NSTEP, rtol=TOL)
    vf = psoln[1]
    return vf
    
#Calculate the final velocity after propagating from the top of the
#atmosphere to the detector, account for all the steps
#Recommend using target="MPI" or "SUF" depending on the detector
def calcVfinal_full(vi, theta,  depth, sigma_p, m_x, interaction="SI", target="full"):
    vf = 1.0*vi
    if (target in ["atmos", "full", "no_shield", "SUF", "MPI", "EDE", "MOD", "surface"]):
        vf = calcVfinal(vf, theta,  depth, sigma_p, m_x, interaction, target="atmos")
    if (target in ["earth", "full", "no_shield", "SUF", "MPI", "EDE", "MOD", "surface"]):
        vf = calcVfinal(vf, theta,  depth, sigma_p, m_x, interaction, target="earth") 
    if (target == "MPI"):
        vf = calcVfinal_shield_MPI(vf, sigma_p, m_x, interaction)
    if (target == "SUF"):
        vf = calcVfinal_shield_SUF(vf, sigma_p, m_x, interaction)
    if (target == "EDE"):
        vf = calcVfinal_shield_EDE(vf, theta, sigma_p, m_x, interaction)
    if (target == "MOD"):
        vf = calcVfinal_shield_MOD(vf, sigma_p, m_x, interaction)
        
    return vf
    
#Calculate the initial velocity (for a given final velocity) after propagating across 'target'
#Here, target = "atmos" or "earth"
def calcVinitial(vf, theta,  depth, sigma_p, m_x, interaction="SI", target="earth"):
    params = [theta, depth, sigma_p, m_x, interaction, target]

    #Propagate across the atmosphere
    if (target == "atmos"):
        d1 = 0
        d2 = pathLength(depth, theta) - pathLength_Earth(depth, theta)

    #Propagate from the surface of the Earth to the detector
    if (target == "earth"):
        d1 = pathLength(depth, theta) - pathLength_Earth(depth, theta)
        d2 = pathLength(depth, theta)
    
    psoln = odeint(dv_by_dD, vf, [d2,d1], args=(params,), mxstep=NSTEP, rtol=TOL)
    return psoln[1]
    
#Calculate the initial speed at the top of the atmosphere for a 
#given final speed at the detector
#Recommend using target="MPI" or "SUF" depending on the detector
def calcVinitial_full(vf, theta,  depth, sigma_p, m_x, interaction="SI", target="full"):
    vi = 1.0*vf
    if (target == "MPI"):
        vi = calcVinitial_shield_MPI(vi, sigma_p, m_x, interaction)
    if (target == "SUF"):
        vi = calcVinitial_shield_SUF(vi, sigma_p, m_x, interaction)
    if (target == "EDE"):
        vi = calcVinitial_shield_EDE(vi, theta, sigma_p, m_x, interaction)
    if (target == "MOD"):
        vi = calcVinitial_shield_MOD(vi, sigma_p, m_x, interaction)
        
    if (target in ["earth", "full", "no_shield", "SUF", "MPI", "EDE", "MOD", "surface"]):
        vi = calcVinitial(vi, theta,  depth, sigma_p, m_x, interaction, target="earth")
    if (target in ["atmos", "full", "no_shield", "SUF", "MPI", "EDE", "MOD", "surface"]):
        vi = calcVinitial(vi, theta,  depth, sigma_p, m_x, interaction, target="atmos")

    return vi
    
#------------------
# Individual function for specifying propagation through shielding at different experiments/sites
    
#Calculate final (or initial) speed after crossing the Lead shielding at SUF
def calcVfinal_shield_SUF(v0, sigma_p, m_x, interaction):
    params = [sigma_p, m_x, interaction]
    #Propagate through 16cm of Lead
    psoln = odeint(dv_by_dD_Pb, v0, [0,16.0e-2] , args=(params,), mxstep=NSTEP, rtol=TOL)
    return psoln[1]
    
def calcVinitial_shield_SUF(v0,  sigma_p, m_x, interaction):
    params = [sigma_p, m_x, interaction]
    #Propagate through 16cm of Lead (backwards)
    psoln = odeint(dv_by_dD_Pb, v0, [16.0e-2,0] , args=(params,), mxstep=NSTEP, rtol=TOL)
    return psoln[1]
    
    
def calcVfinal_shield_MOD(v0, sigma_p, m_x, interaction):
    params = [sigma_p, m_x, interaction]
    #Propagate through 20cm of Lead
    psoln = odeint(dv_by_dD_Pb, v0, [0,20.0e-2] , args=(params,), mxstep=NSTEP, rtol=TOL)
    return psoln[1]
    
def calcVinitial_shield_MOD(v0,  sigma_p, m_x, interaction):
    params = [sigma_p, m_x, interaction]
    #Propagate through 20cm of Lead (backwards)
    psoln = odeint(dv_by_dD_Pb, v0, [20.0e-2,0] , args=(params,), mxstep=NSTEP, rtol=TOL)
    return psoln[1]
    


#Calculate final (or initial) speed after crossing the Copper shielding at MPI
def calcVfinal_shield_MPI(v0, sigma_p, m_x, interaction):
    params = [sigma_p, m_x, interaction]
    #Propagate through 1mm Copper
    psoln = odeint(dv_by_dD_Cu, v0, [0,1e-3] , args=(params,), mxstep=NSTEP, rtol=TOL)
    return psoln[1]
    
def calcVinitial_shield_MPI(v0, sigma_p, m_x, interaction):
    params = [sigma_p, m_x, interaction]
    #Propagate through 1mm Copper
    psoln = odeint(dv_by_dD_Cu, v0, [1e-3,0] , args=(params,), mxstep=NSTEP, rtol=TOL)
    return psoln[1]

hdet = 0.3 #height of EDE detector from the inner lead radius, 30cm
hroom = 3.4 #height of room, 3.4m
hfloor = 1.5 #height of EDE detector from the floor of the room, 1.5m
Lroom = 1.5 #Side length of the room, 4m 

a_ceil = np.pi - np.arctan(Lroom/(hroom - hfloor))
a_floor = np.arctan(Lroom/hfloor)
a_bottom = np.arctan(0.5/(hdet + 0.1))
a_lead = np.pi - np.arctan(0.234/(0.83 - hdet))
a_steel = np.pi - np.arctan(0.215/(1.5-hdet))

#Calculate final (or initial) speed after crossing the Lead shielding at EDE...
def calcVfinal_shield_EDE(v0, theta, sigma_p, m_x, interaction):
    params = [sigma_p, m_x, interaction]

    #Walls of the room                                                                                                                                       
    D_walls = 0.4
    if (theta > a_ceil):
        D_walls = 0.2
    elif (theta <a_floor):
        D_walls= 0.1

    delta = 0.0
    D_walls *= (1+delta)

    v1 = odeint(dv_by_dD_concrete, v0, [0,D_walls] , args=(params,), mxstep=NSTEP, rtol=TOL)[1]

    #Steel bottom
    if (theta < a_bottom):
        v1 = odeint(dv_by_dD_Fe, v1, [0,10.e-2] , args=(params,), mxstep=NSTEP, rtol=TOL)[1]

    #Lead shielding
    if (theta < a_lead):
        v1 = odeint(dv_by_dD_Pb, v1, [0,10.e-2] , args=(params,), mxstep=NSTEP, rtol=TOL)[1]

    #Steel can
    D_steel = 1.5e-3
    if (theta > a_steel):
        D_steel = 2.5e-2
    v1 = odeint(dv_by_dD_Fe, v1, [0,D_steel] , args=(params,), mxstep=NSTEP, rtol=TOL)[1]

    #Copper elements
    D_copper = 6e-2
    v1 = odeint(dv_by_dD_Cu, v1, [0,D_copper] , args=(params,), mxstep=NSTEP, rtol=TOL)[1]

    return v1
    
def calcVinitial_shield_EDE(v0, theta,  sigma_p, m_x, interaction):
    params = [sigma_p, m_x, interaction]
 
    #Copper elements                                                                                                                                         
    D_copper = 6e-2
    v1 = odeint(dv_by_dD_Cu, v0, [D_copper,0] , args=(params,), mxstep=NSTEP, rtol=TOL)[1]

    #Steel can                                                                                                                                               
    D_steel = 1.5e-3
    if (theta > a_steel):
        D_steel = 2.5e-2
    v1 = odeint(dv_by_dD_Fe, v1, [D_steel,0] , args=(params,), mxstep=NSTEP, rtol=TOL)[1]

    #Lead shielding                                                                                                                                          
    if (theta < a_lead):
        v1 = odeint(dv_by_dD_Pb, v1, [10.e-2,0] , args=(params,), mxstep=NSTEP, rtol=TOL)[1]

    #Steel bottom                                                                                                                                            
    if (theta < a_bottom):
        v1 = odeint(dv_by_dD_Fe, v1, [10.e-2,0] , args=(params,), mxstep=NSTEP, rtol=TOL)[1]

    #Walls of the room                                                                                                                                       
    D_walls = 0.4
    if (theta > a_ceil):
        D_walls = 0.2
    elif (theta < a_floor):
        D_walls = 0.1

    v1 = odeint(dv_by_dD_concrete, v1, [D_walls,0] , args=(params,), mxstep=NSTEP, rtol=TOL)[1]

    return v1

    
