import numpy as np
from numpy import pi
from scipy.integrate import simps,quad, cumtrapz
import verne
from LabFuncs import *
import utils
from scipy.special import erf
from scipy.interpolate import interp1d
import MaxwellBoltzmann as MB
import argparse
import os.path

import CalcVelDist_function

import matplotlib.pyplot as plt

try:
  import WIMpy
  WIMpy_installed = True
  
except ImportError:
  print(">----- WARNING: WIMpy is not installed. -----<")
  print(">----- To include checks using WIMpy, install at https://github.com/bradkav/WIMpy_NREFT -----<")
  WIMpy_installed = False


#Reference values
#sigma_ref = 1e-32
#m_ref = 0.5

#Number of grid points in gamma and v
N_gamma = 11
N_v = 62
gamma_grid = np.linspace(0, 1, N_gamma)*np.pi

#Load velocity distribution from file
def getVelDist(m_x, sigma, gamma_ind):
    
    gamma, v, f = np.loadtxt("../results/veldists/f_SI_surface_mx%.3f_lsig%.2f.txt"%(m_x, np.log10(sigma)), unpack=True)
    
    #vvals = vvals1[gamma_ind*61:(gamma_ind+1)*61]
    #fvals = fvals1[gamma_ind*61:(gamma_ind+1)*61]
    
    v_grid = np.reshape(v, (N_gamma, N_v))
    f_grid = np.reshape(f, (N_gamma, N_v))
    #gamma_list = gamma[::62]
    
    return v_grid[gamma_ind,:], f_grid[gamma_ind,:]

#Calculate velocity integral from an interpolation function
#defining f(v) and a maximum speed vmax=800 km/s
def calcEta(v, interpfun):
    return quad(lambda x: interpfun(x)/x, v, 800.0, epsrel=1e-3)[0]

#Prefactor with correct normalisation for DD signal
def rate_prefactor(A, m_x):
    rho0 = 0.3
    mu = 1.78e-27*(m_x*0.9315)/(m_x + 0.9315)
    
    #Here, there was an extra factor of 1/(2pi) which shouldn't be there.
    #return 1.38413e-12*rho0/(4.0*np.pi*m_x*mu*mu)
    return 1.38413e-12*rho0/(2.0*m_x*mu*mu)

#Calculate recoil spectrum
@np.vectorize
def dRdE(E, A, mx,sig,f_interp):  
    int_factor = sig*verne.calcSIFormFactor(E, A)*A**2
    return rate_prefactor(A, mx)*int_factor*calcEta(MB.vmin(E, A, mx),  f_interp)


#-------------------------------------------------------------

#Filename of test velocity distribution
fv_file = "../results/veldists/f_SI_surface_mx0.500_lsig-36.00.txt"

#Check whether the test velocity distribution exists, 
#otherwise, generate it

if  (os.path.isfile(fv_file) != True):
    print(">Test velocity distribution <../results/veldists/f_SI_surface_mx0.500_lsig-36.00.txt> not found.")
    print(">Generating from scratch...")
    #os.system("python3 CalcVelDist.py -m_x 0.5 -sigma_p 1e-36 -loc surface -int SI")
    CalcVelDist_function.calcVelDist_full(m_x=0.5, sigma_p=1e-36, loc="surface", interaction = "SI")
    print(">...done.")
else:
    print(">Loading test velocity distribution <../results/veldists/f_SI_surface_mx0.500_lsig-36.00.txt>.")



print(">Checking normalisation, int f(v) dv")
print("gamma/pi \t normalisation")
print("-------- \t -------------")


for i in range(N_gamma):
    #Define an interpolation function for f(v)
    vvals,fvals = getVelDist(0.5, 1e-36, gamma_ind=i)
    #Check that the distribution function is properly normalised
    norm = np.trapz(fvals, vvals)
    print("%.2f \t\t %.4f"%(gamma_grid[i]/np.pi, norm))
    
    
#Check the rates for gamma = pi
vvals,fvals = getVelDist(0.5, 1e-36, gamma_ind=(N_gamma-1))
f_interp = interp1d(vvals, fvals, kind='linear',bounds_error=False, fill_value=0.0)
    

E_list = np.geomspace(1e-3, 120e-3, 1000) #in keV
dRdE_verne = dRdE(E_list, 28, 0.500, 1e-36, f_interp)


#If WIMpy is installed - cross check using that code too
if (WIMpy_installed):
    #Check against the standard rate without Earth-stopping
    dRdE_WIMpy_free = WIMpy.DMUtils.dRdE_standard(E_list, N_p=14, N_n=14, m_x=0.5, sig=1e-36, vlag=232, vesc=544.0)
    
    #WIMpy also allows you to feed in a generic function for
    # eta(vmin) = int_vmin^infty f(v)/v dv, so let's try that
    eta_list = np.trapz(fvals/vvals, vvals) - cumtrapz(fvals/vvals, vvals, initial=0)
    eta_interp = interp1d(vvals, eta_list, bounds_error=False, fill_value = 0.0)
    
    dRdE_WIMpy_verne = WIMpy.DMUtils.dRdE_generic(E_list, N_p=14, N_n=14, m_x=0.5, sig=1e-36, eta_func=eta_interp)
    

plt.figure()

plt.semilogy(E_list*1e3, dRdE_verne, label="Verne", lw=2)

if (WIMpy_installed):
    plt.semilogy(E_list*1e3, dRdE_WIMpy_free, label="WIMpy (no Earth-scattering)", lw=2)
    plt.semilogy(E_list*1e3, dRdE_WIMpy_verne, color='k', linestyle=':', label="WIMpy (verne velocity distribution)", lw=2)

plt.xlim(0, 120)
plt.ylim(1e1, 1e7)

plt.xlabel("$E_R$ [eV]")
plt.ylabel("$\mathrm{d}R/\mathrm{d}E_R$ [1/keV/kg/day]")

plt.title(r"$m_\chi = 500 \,\mathrm{MeV}; \sigma_p^{\mathrm{SI}} = 10^{-36}\,\mathrm{cm}^2$; $\gamma = \pi$; Silicon target")

plt.legend()

plt.savefig("../plots/TestSpectrum.pdf", bbox_inches='tight')
plt.show()

