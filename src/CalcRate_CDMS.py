import numpy as np
from numpy import pi
from scipy.integrate import simps,quad
import verne
from LabFuncs import *
import utils
from scipy.special import erf
from scipy.interpolate import interp1d
import MaxwellBoltzmann as MB
import argparse
import os.path


#Parse the arguments!
parser = argparse.ArgumentParser(description='...')
parser.add_argument('-m_x','--m_x', help='DM mass in GeV', type=float,default = 1e5)
parser.add_argument('-sigma_p','--sigma_p', help='DM-nucleon cross section, sigma_p in cm^2', type=float, required=True)
args = parser.parse_args()
m_x = args.m_x
sigma_p = args.sigma_p


#guesstimating from Fig. 8 of CDMS-I paper
#FWHM energy resolution at threshold is about 1 keV
FWHM = 1
sig_th = 2.355*FWHM

#Exposure in kg days
exposure = 15.8

#CDMS-I efficiency function
def eff(E):
    E1 = 10.0
    E2 = 100.0
    return 0.5*(erf((E-E1)/(np.sqrt(sig_th))) + erf((E2-E)/(np.sqrt(sig_th))))

#Prefactor with correct normalisation for DD signal
def rate_prefactor(A, m_x):
    rho0 = 0.3
    mu = 1.78e-27*(m_x*0.9315)/(m_x + 0.9315)
    return 1.38413e-12*rho0/(2.0*m_x*mu*mu)

#Stanford
lat_SUF = +37.4 #N
lon_SUF = -122.2 #W

#Get Julian date of exposure
t0 = JulianDay(11, 1, 1998, 1)
Nvals = 10001

#Pick a typical day to integrate over
#April 21st...
t1 = 171.0
t2 = 172.0
tvals = t0 + np.linspace(t1, t2, Nvals)

gammavals = np.zeros(Nvals)

#Calculate gamma from the LabVelocity
for i in range(Nvals):
    vs = -LabVelocity(tvals[i], lat_SUF, lon_SUF)
    vs_hat = vs/np.sqrt(np.sum(vs**2))
    rdet_hat = np.asarray([0,0,1.0])
    gammavals[i] = np.arccos(np.dot(vs_hat, rdet_hat))

#Load velocity distribution from file
def getVelDist(mstr, lsigstr, gamma_ind):
    Ngamvals = 11
    Nvvals = 61

    gamma_vals1, vvals1, fvals1 = np.loadtxt("../results/veldists/f_SUF_lmx" + mstr+"_lsig" + lsigstr + ".txt", unpack=True)
    vvals = vvals1[gamma_ind*61:(gamma_ind+1)*61]
    fvals = fvals1[gamma_ind*61:(gamma_ind+1)*61]
    return vvals, fvals

#Calculate velocity integral from an interpolation function
#defining f(v) and a maximum speed vmax
def calcEta_final(v, interpfun, vmax):
    return quad(lambda x: interpfun(x)/x, v, vmax*1.1)[0]

#Calculate recoil spectrum
def dRdE(E, A, mx,sig,interpfun,vmax):  
    int_factor = sig*verne.calcSIFormFactor(E, A)*A**2
    return rate_prefactor(A, mx)*int_factor*calcEta_final(MB.vmin(E, A, mx),  interpfun, vmax)

#Calculate number of signal events
def Nevents(E_min, E_max, m_x, sig, gamma_ind = 10):
    sigstring = '{0:.2f}'.format(np.log10(sig))
    mstring = '{0:.1f}'.format(np.log10(m_x))
    
    vvals,fvals = getVelDist(mstring, sigstring, gamma_ind)
    vmax = np.max(vvals)

    interpfun = interp1d(vvals, fvals, kind='linear',bounds_error=False, fill_value=0.0)
    integ = lambda x: eff(x)*dRdE(x, 73.0, m_x, sig, interpfun,vmax)
    return quad(integ, E_min, E_max, epsrel=1e-4)[0]



#Calculate number of events as a function of gamma
Ne_list = np.zeros(11)
print "gamma    Ne"
for i in range(11):
    #We integrate from 1 -> 110 keV because the efficiency function takes care of the thresholds...
    Ne_list[i] = Nevents(1.0, 110.0, m_x, sigma_p, gamma_ind=i)*exposure
    print i*np.pi/10.0, Ne_list[i]

Ne_interp = interp1d(np.linspace(0, 1,11)*np.pi, Ne_list)

#Integrate over the values of gamma for a single day
Ne_tot =  simps(Ne_interp(gammavals), tvals-t0)
print "Total number of events:", Ne_tot

#Append to number of events file
fname = "../results/Nevents/N_SUF_lmx" + '{0:.1f}'.format(np.log10(m_x)) + ".txt"
outarray = np.c_[np.log10(sigma_p), Ne_tot]

if (not os.path.isfile(fname)):
    htxt = "Total number of CDMS-I signal events at SUF (signal averaged over 1 day). Log10(m_x) = "+"{0:.1f}".format(np.log10(m_x))+"\nColumns: Log10(sigma/cm^2)      N_sig"
    np.savetxt(fname, outarray, header=htxt)
else:
    f_handle = file(fname, 'a')
    np.savetxt(f_handle, outarray)
    f_handle.close()



