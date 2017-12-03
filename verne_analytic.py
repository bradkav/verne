import numpy as np
import WIMpy.DMUtils as DMU
from numpy import pi
from scipy.integrate import quad, dblquad, tplquad, fixed_quad, quadrature, romberg
from scipy.interpolate import interp1d, interp2d, griddata, RectBivariateSpline,InterpolatedUnivariateSpline
import os
import sys
import scipy.special
from scipy.special import erf
from skmonaco import mcquad
from scipy.integrate import odeint

from tqdm import tqdm

import mcint

import matplotlib.pyplot as pl

ROOT2 = np.sqrt(2.0)

vcut = 0.1

#15cm of Lead + 1cm of ancient lead (A = 207)
#Number density 3.3e22 per cc

#about 3cm copper A = 63.5
#Number density 8.5e22 per cc

#25cm of polyethylene
#Number density 

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
Niso_full = None
Deff_interp = None
phi_interp = None
corr_interp = None
xmax = np.zeros(8)
corr_Pb = None


Hvals = None
Tvals = None
beta = None

isoID = {"O":0, "Si":1, "Mg":2, "Fe":3, "Ca":4, "Na":5, "S":6, "Al":7, "O_A":8, "N_A": 9}

#h_A = 80e3
h_A = 1e-5
print " - CURRENTLY NO ATMOSPHERE!!!"
print " - CURRENTLY NO ATMOSPHERE!!!"
print " - CURRENTLY NO ATMOSPHERE!!!"
print " - CURRENTLY NO ATMOSPHERE!!!"
print " - CURRENTLY NO ATMOSPHERE!!!"
print " - CURRENTLY NO ATMOSPHERE!!!"

print "Carefully separate shield, atmosphere and Earth!"

R_E = 6371.0e3 + h_A
#R_A = R_E + 1e5

#Densities are in number/cm^3
#Distances in m



def loadIsotopes():
    print "    Loading isotope data and density profiles..."
    
    global dens_profiles
    global isotopes
    global Avals
    global dens_interp
    global Niso
    global Niso_full
    global corr_interp
    global corr_Pb
    
    global Hvals
    global Tvals
    global beta
    
    rootdir = "data/"
    
    Avals = np.loadtxt(rootdir+"isotopes.txt", usecols=(1,))[:-1] #Trim off Xenon
    isotopes = np.loadtxt(rootdir+"isotopes.txt", usecols=(0,))[:-1] #Trim off Xenon
    Niso = len(isotopes)
    Niso_full = Niso + 2
    #print Niso_full
    
    r_list = np.loadtxt(rootdir+"dens_profiles/n_1.dat", usecols=(0,))
    r_list0 = 0.0*r_list
    dens_profiles = [np.loadtxt(rootdir+"dens_profiles/n_"+str(int(iso))+".dat", usecols=(1,)) for iso in isotopes]
    
    h_list = np.linspace(0, h_A, 100)+1e-5
    
    #Make a slight correction of the profile - truncate at 6371 km
    r_list[-1] = R_E - h_A
    r_list = np.append(r_list, R_E - h_A + h_list)
    #print r_list
    #print Niso_full
    for i, dens in enumerate(dens_profiles):
        dens[-1] = dens[-2]
        dens_profiles[i] = np.append(dens, 0.0*h_list)
    
    #Add the atmospheric elements:
    Avals = np.append(Avals,[16, 14])
    #print Avals
    #Calculate the atmosphere density profiles...
    Hvals, Tvals, beta = np.loadtxt(rootdir+"ISA/temp.txt", unpack=True)
    Hvals *= 1e3
    beta *= 1e-3 #Get everything in m
    
    frac = [0.21, 0.78]
    dens_vec = np.vectorize(density)
    dens_atmos = [np.append(r_list0, 2*f_n*dens_vec(h_list)) for f_n in frac]
    dens_profiles.extend(dens_atmos)
    
    
    dens_interp = [interp1d(r_list, dens, bounds_error=False, fill_value = 0) for dens in dens_profiles]
    #print dens_interp[0]
    corr_interp = [calcFFcorrection(Avals[ID]) for ID in range(Niso_full)]
    
    print " Need to carefully check the density profiles..."
    corr_Pb = calcFFcorrection(207)
    
    """
    pl.figure()
    for i in range(Niso_full):
        pl.plot(r_list, dens_interp[i](r_list), '-',label=str(i))
    pl.legend(loc='best')
    pl.show()
    """
    
def density(h):
    H = (R_E-h_A)*h/((R_E-h_A) + h)
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
    
def pathLength(depth, theta):
    #Theta = 0 is directly from BELOW, angles in radians
    r_det = R_E - depth - h_A
    return +np.cos(theta)*r_det + np.sqrt((-np.cos(theta)*r_det)**2 - (r_det**2 - R_E**2))
    
def calcDeff(depth,theta,ID):
    r_det = R_E - depth - h_A
    d_total = pathLength(depth, theta)
    res = quad(lambda l: dens_interp[ID](np.sqrt(r_det**2 + l**2 - 2*l*r_det*np.cos(theta))), 0, d_total)[0]
    #Convert from m*cm^-3 to cm^-2
    return 1e2*res
    
def calcDeff_interp(depth, ID):
    global Deff_interp
    
    tvals = np.linspace(0, np.pi,100)
    Deff_vals = tvals*0.0
    for i, t in enumerate(tvals):
        Deff_vals[i] = calcDeff(depth,t,ID)
    Deff_interp = interp1d(tvals, Deff_vals, kind="linear")

def calcFFcorrection(A0):
    v_vals = np.linspace(0, 1000, 200)

    #print v_vals
    corr_fact = v_vals*0.0
    for i, v in enumerate(v_vals):
            corr_fact[i] = quad(lambda x: 2.0*x*DMU.calcSIFormFactor(x*(1e6/(3e5*3e5))*ERmax(1e10, 0.9315*A0, v), A0), 0, 1)[0]
    corr_fact[0] = 1.0
    return interp1d(v_vals, corr_fact, kind='linear', bounds_error=False, fill_value=0.0)



def effectiveXS(sigma_p, m_X, ID):
    return sigma_p*(0.9315/m_X)*Avals[ID]**5
    
    
    
def CalcFullF_full(v, gamma, depth,sigma_p, m_x, target, tcut):
    integrand = lambda t: f_integrand_full(v, t, gamma, depth, sigma_p, m_x, target)
    
    #Calculate cut off velocity    
    #return romberg(integrand, tcut*1.0, np.pi, rtol=1e-2, vec_func=False,show=1)
    #return quadrature(integrand, tcut, np.pi, tol=1e-1,rtol=1e-4)[0] #THIS ONE WORKS OKAY!
    #return fixed_quad(integrand, tcut, np.pi)[0] 
    #epsabs=1e-3, epsrel=1e-3
    integ_res = quad(integrand, tcut, np.pi, epsabs=0, epsrel=1e-1)
    print "Result:", integ_res[0], "+-", integ_res[1]
    #print " "
    return integ_res[0] #This seems to work too!
    #Be careful about what range of values we care about here for theta
    
def CalcFullEta(vmin, gamma, sigma_p, m_x, ID, vb=760.0):
    #print calcVfinal(800.0, np.pi/2.0, sigma_p, m_x, ID)
    integrand = lambda x: eta_integrand(x[0], x[1], gamma, sigma_p, m_x, ID)
    res, error = mcquad(integrand, npoints=1e4, xl=[vmin,0.0], xu=[vb,np.pi])
    return res, error
    

def f_integrand_full(vf, t1, gamma, depth, sigma_p, m_x, target="earth"):
    #if (cut(theta) < 1.0):
    #    return 0.0
    #print "HERE!"
    #vi = calcVinitial(vf, theta, sigma_p, m_x,  ID)
    if hasattr(t1, "__len__"):
        theta = t1[0]
    else: 
        theta = t1
    dv = 1e-1

    #try:
    vi1 = calcVinitial_full(vf+dv/2.0, theta,  depth, sigma_p, m_x, target)
    vi2 = calcVinitial_full(vf-dv/2.0, theta,  depth, sigma_p, m_x, target)
       #print vi1, vi2
    #except ValueError as err:
       #print err.args
    #   print "Exiting..."
       #print "High v:", vf
    #   return 0.0
    #print sigma_p, vi, vf
    #print " This isn't right!"
    #Estimate the derivative
    
    vi = (vi1 + vi2)/2.0
    deriv = np.abs(vi1 - vi2)*1.0/dv
    
    if (vi > 790):
        return 0
    
    #deriv = np.abs(calcVinitial(vf+dv/2.0, theta, sigma_p, m_x,  ID) - calcVinitial(vf-dv/2.0, theta, sigma_p, m_x,  ID))*1.0/dv
    #print deriv
    #deriv = 1.0
    
    #deriv = vi/vf

    
    #if (deriv < 1e-10):
    #    return 0.0
        
    res = (deriv)*np.sin(theta)*(vi**2)*calcf_integ(vi, theta, gamma)
    #print vf, theta, res
    #print vf, theta/np.pi,vi1,  res
    #In principle, should be vf**2/deriv...
    return res
 
 
def radius(D, theta, depth):
    r_det = R_E - depth - h_A
    return np.sqrt(R_E**2 + D**2 + 2*D*(r_det*np.cos(theta) - pathLength(depth, theta)))
    #return depth + np.sqrt(r_det**2 + D**2 + 2*D*r_det*np.cos(theta))

def f(v, D, params):
    #print v, D
    #print v, (D - h_A)
    #if (v > 800.0):
    #    raise ValueError('A very specific bad thing happened')
    #print v
    #if (v < vcut):
    #   raise ValueError('Cut-off value reached')
    theta, depth, sigma_p, m_x, target = params
    res = 0.0
    if (target == "atmos"):
        isovals = [8,9]
    elif (target == "earth"):
        isovals = range(Niso)
    else:
        isovals = range(Niso_full)
    
    r = radius(D, theta, depth)
    #if (r > R_E):
    #    print "Help!", r, R_E
    for i in isovals:
        #print i
        res += dens_interp[i](r)*effectiveXS(sigma_p, m_x, i)*corr_interp[i](v)
    return -1e2*v*res #(km/s)/m

    
def f_back(v, D, params):
    return -f(v,D,params)
    
def f2(v, D, params):
    if (v < 0):
        return 0
        
    #Pb density
    n_Pb = 3.3e22
    A0 = 207
        
    theta, depth, sigma_p, m_x = params
    res = n_Pb*sigma_p*(0.9315/m_x)*A0**5*corr_Pb(v)
    return -1e2*v*res #(km/s)/m

def calcVfinal_full(v0, theta,  depth, sigma_p, m_x, target="full"):
    params = [theta, depth, sigma_p, m_x, target]
    #print pathLength(depth, theta)
    

    d = depth
    if (target == "atmos"):
        d = 0
        
    #print theta, pathLength(d, theta)
        
    try:
        psoln = odeint(f, v0, [0,pathLength(d, theta)] , args=(params,), mxstep=1000, rtol=1e-6)
    except ValueError as err:
        #print v0, theta
        return vcut
    #print v0, theta
    return psoln[1]
    #Initial values
    
def calcVinitial_full(v0, theta,  depth, sigma_p, m_x, target="earth"):
    params = [theta, depth, sigma_p, m_x, target]
    #PL = pathLength(depth, theta)
    #print PL
    #lam = sigma_p*1e35
    #print 1e4*(m_x/lam)
    #if (PL - h_A > 1e4*(m_x/lam)):
    #    return 799.0
    #Need to carefully check/think about this criterion
    
    
    
    #Calculate up to surface
    d = depth
    if (target == "atmos"):
        d = 0
        
    #print theta, pathLength(d, theta)
    #print "Line 429: Make sure this is symmetric with the calcVfinal_full function..."
    
    #psoln = odeint(f, v0, [pathLength(depth, theta),pathLength(0, theta)], args=(params,), mxstep=1000)
    psoln = odeint(f, v0, [pathLength(d, theta),0], args=(params,), mxstep=1000, rtol=1e-4)
    #print v0, psoln[1]
    return psoln[1]
    #Initial values
    
def calcVfinal_shield(v0, theta,  depth, sigma_p, m_x):
    params = [theta, depth, sigma_p, m_x]
    #print pathLength(depth, theta)
    psoln = odeint(f2, v0, [0,16.0e-2] , args=(params,), mxstep=1000)
    return psoln[1]
    
print " I should be able to calculate the atmosphere separately! - Double-check everything"
    
    
#------- Velocity distribution stuff----------
#----------------------------------------------

def IntegralOverPhi(x, phi_max):
    if (phi_max < 0):
        return 0
    if (phi_max >= np.pi):
        return 2.0*np.pi*scipy.special.i0(x)
    else: 
        return 2.0*phi_interp(x, phi_max)

IntegralOverPhiVec = np.vectorize(IntegralOverPhi)

def calcf_integ(v, theta, gamma, sigmav=156.0, vesc=533.0):
    
    if (np.sin(gamma)*np.sin(theta) <= 1e-10):
        return 2.0*np.pi*VelDist(v, theta, 0, gamma)
    
    ve = ROOT2*sigmav
    delsq = v**2 + ve**2 - 2*v*ve*np.cos(gamma)*np.cos(theta)
    
    cosmin = (v**2 + ve**2 - vesc**2)/(2*v*ve*np.sin(gamma)*np.sin(theta))\
         - (np.cos(gamma)*np.cos(theta))/(np.sin(gamma)*np.sin(theta))
    
    
    x0 = np.sin(theta)*np.sin(gamma)*v*ve/(sigmav**2)
    phi_max = np.arccos(np.clip(cosmin, -1.0, 1.0))
    #print cosmin, phi_max/np.pi
    A = IntegralOverPhiVec(x0, phi_max)*np.exp(-delsq/(2.0*sigmav**2))
    
    return A*1.0/NNORM
    
def VelDist(v, theta, phi, gamma, sigmav=156.0, vesc=533.0):
    cdel = np.sin(gamma)*np.sin(theta)*np.cos(phi) + np.cos(gamma)*np.cos(theta)
    ve = ROOT2*sigmav
    dsq = v**2 - 2*v*ve*cdel + ve**2
    A = np.exp(-dsq/(2*sigmav**2))/NNORM
    if hasattr(A, "__len__"):
        A[np.where(dsq > vesc**2)] = A[np.where(dsq > vesc**2)]*0.0
    else:
        if (dsq > vesc**2):
            A = 0
    return A
    
# N_esc - normalisation constant
def Nesc(sigmav=156.0, vesc=533.0):
    return (erf(vesc/(ROOT2*sigmav)) - np.sqrt(2.0/pi)*(vesc/sigmav)*np.exp(-vesc**2/(2.0*sigmav**2)))

NNORM = Nesc()*156.0**3*np.sqrt(2.0*np.pi)*2.0*np.pi
    
def loadPhiInterp():
    global phi_interp
    fname = "data/PhiIntegrals.dat"
    x = np.arange(-7, 7.001, 0.05)
    y = np.arange(0, np.pi+0.1, 0.05)
    xlen = len(x)
    ylen = len(y)
    #print xlen, ylen
    if (os.path.isfile(fname)):
        data = np.loadtxt(fname, usecols=(2,))
        z = data.reshape((xlen,ylen))
        phi_interp = interp2d(x,y,z.T)

    else:
        print " File 'data/PhiIntegrals.dat' doesn't exist..."
        print " Exiting..."
        sys.exit()
    