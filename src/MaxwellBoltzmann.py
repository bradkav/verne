import numpy as np
from scipy.integrate import quad
try:
    from scipy.integrate import simps
except ImportError:
    from scipy.integrate import simpson as simps
    
from scipy.interpolate import interp2d, RectBivariateSpline
import scipy.special
import os

from tqdm import tqdm
    
phi_interp = None

#------- Velocity distribution stuff----------
#----------------------------------------------


vesc = 544.0 #Escape velocity
sigmav = 156.0#238/np.sqrt(2) #156.0  #Velocity dispersion (v_0/sqrt(2))
ve = 263.0#263.0 #232.0  #Earth peculiar velocity around GC


# Nesc - normalisation constant
Nesc = (scipy.special.erf(vesc/(np.sqrt(2.0)*sigmav)) - np.sqrt(2.0/np.pi)*(vesc/sigmav)*np.exp(-vesc**2/(2.0*sigmav**2)))
#NNORM = Nesc*156.0**3*np.sqrt(2.0*np.pi)*2.0*np.pi
NNORM = Nesc*sigmav**3*np.sqrt(2.0*np.pi)*2.0*np.pi



#Load an interpolation function for the integral
#of the Maxwell-Boltzmann velocity distribution
#over phi
#*This is called as soon as the module is loaded*
def loadPhiInterp():
    global phi_interp
    curr_dir = os.path.dirname(os.path.realpath(__file__)) + '/'
    fname = curr_dir + "../data/PhiIntegrals.dat"
    xvals = np.linspace(-10, 10, 1001)
    phivals = np.linspace(0, np.pi, 501)
    xlen = len(xvals)
    ylen = len(phivals)

    if (os.path.isfile(fname)):
        data = np.loadtxt(fname, usecols=(2,))
        z = data.reshape((xlen,ylen))

    else:
        print(">MAXWELLBOLTZMANN: File '../data/PhiIntegrals.dat' doesn't exist...")
        print(">MAXWELLBOLTZMANN: Calculating from scratch...")
        z = np.zeros((xlen, ylen))
        for i,x in enumerate(tqdm(xvals)):
            for j,phi in enumerate(phivals):
                #if (j == 0):
                #    continue
                #_phis = np.linspace(0, phi, 101)
                #z[i,j] = simps(np.exp(x*np.cos(_phis)), _phis)
                z[i,j] = quad(lambda y: np.exp(x*np.cos(y)), 0, phi, epsabs=1e-12)[0]
                #z[i,j] = np.clip(z[i,j], None, np.pi*scipy.special.i0(x))

        xgrid, phigrid = np.meshgrid(xvals, phivals, indexing='ij')
        np.savetxt(fname, np.c_[xgrid.flatten(), phigrid.flatten(),z.flatten()])
    
    phi_interp = RectBivariateSpline(xvals, phivals, z, kx = 1, ky = 1)

loadPhiInterp()

def IntegralOverPhi(x, phi_max):
    if (np.abs(x) > 10):
        print(x, phi_max)
        assert 1 == 0
    if (phi_max <= 0):
        return 0
    if (phi_max >= np.pi):
        return 2.0*np.pi*scipy.special.i0(x)
    else: 
        return 2.0*phi_interp(x, phi_max)
    

def IntegralOverPhiVec(x, phi_max):
    if hasattr(phi_max, "__len__"):
        if (np.sum(np.abs(x) > 10) > 0):
            raise ValueError("Values of x in IntegralOverPhiVec exceeds tabulated values. Retabulate!")
        result = np.zeros(len(phi_max))
        inds = phi_max >= np.pi
        result[inds] = 2.0*np.pi*scipy.special.i0(x[inds])

        inds = (0 <= phi_max)  & (phi_max <= np.pi)
        result[inds] = 2.0*phi_interp(x[inds], phi_max[inds], grid=False)
        return result
    else:
        return IntegralOverPhi(x, phi_max)
    
    
#Integrand for integrating over the velocity distribution
#The phi integral has already been performed, so all you have
#left is v and theta.
def calcf_integ(v, theta, gamma):
    #print("Vectorizing!")
    
    fudge = 1e-6
    theta = np.clip(theta, fudge*np.pi, (1-fudge)*np.pi)
    gamma = np.clip(gamma, fudge*np.pi, (1-fudge)*np.pi)
    
    sin_theta = np.sin(theta) #+ 1e-6
    sin_gamma = np.sin(gamma)
    
    #This function is vectorized over `theta` but not over `v`
    if (hasattr(v, "__len__")):
        raise ValueError("Velocity v must be a scalar.")
            
    if (v < 1):
        return 2.0*np.pi*VelDist(v, theta, 0, gamma)
        
    if (v > (ve + vesc)):
        return 0.0
            
    delsq = v**2 + ve**2 - 2*v*ve*np.cos(gamma)*np.cos(theta)
                    
    x0 = sin_theta*sin_gamma*v*ve/(sigmav**2)

    C1 = (v**2 + ve**2 - vesc**2)/(2*v*ve)
    C2 = (np.cos(gamma)*np.cos(theta))

    cosmin = (C1 - C2)/(sin_theta*sin_gamma)

    phi_max = np.arccos(np.clip(cosmin, -1.0, 1.0))

    #print(x0, phi_max/np.pi)
    A = IntegralOverPhiVec(x0, phi_max)*np.exp(-delsq/(2.0*sigmav**2))
    
    return A*1.0/NNORM
    
#Full 3-D velocity distribution (v, theta, phi)
def VelDist(v, theta, phi, gamma):
    cdel = np.sin(gamma)*np.sin(theta)*np.cos(phi) + np.cos(gamma)*np.cos(theta)
    dsq = v**2 - 2*v*ve*cdel + ve**2
    A = np.exp(-dsq/(2.0*sigmav**2))/NNORM
    if hasattr(A, "__len__"):
        A[np.where(dsq > vesc**2)] = A[np.where(dsq > vesc**2)]*0.0
    else:
        if (dsq > vesc**2):
            A = 0
    return A

#Calculate the free MB speed distribution 
#after integrating over all angles
def calcf_SHM(v):
    beta = ve/(sigmav**2)
    N1 = 1.0/(Nesc*sigmav**3*np.sqrt(2*np.pi))
    f = v*0.0
    
    a = (v <= vesc-ve)
    f[a] = np.exp(-(v[a]**2 + ve**2)/(2.0*sigmav**2))*(np.exp(beta*v[a])-np.exp(-beta*v[a]))
    
    b = (vesc-ve < v) & (v < vesc+ve)
    f[b] = np.exp(-(v[b]**2 + ve**2)/(2.0*sigmav**2))*(np.exp(beta*v[b]) - np.exp((v[b]**2 + ve**2 -vesc**2)/(2*sigmav**2)))
    return f*v*N1/beta

#Minimum velocity required for a recoil of energy E_R
def vmin(E, m_N, m_x):
    res = E*0.0
    m_N2 = m_N*0.9315
    mu = (m_N2*m_x)/(m_N2+m_x)
    res =  3e5*np.sqrt((E/1e6)*(m_N2)/(2*mu*mu))
    return res
