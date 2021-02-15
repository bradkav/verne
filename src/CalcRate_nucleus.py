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

#Exposure in kg days
exposure = 0.046e-3

#nucleus efficiency function
def eff(E):
    E_th = 19.7e-3
    sig_th = 3.83e-3
    return 0.5*(1+erf((E - E_th)/(np.sqrt(2)*sig_th)))

#Prefactor with correct normalisation for DD signal
def rate_prefactor(A, m_x):
    rho0 = 0.3
    mu = 1.78e-27*(m_x*0.9315)/(m_x + 0.9315)
    return 1.38413e-12*rho0/(2.0*m_x*mu*mu)

#MPI Munich
lat_MPI = +48.1 #N
lon_MPI = +11.57 #W




#From https://arxiv.org/src/1707.06749v4/anc/additional_material.txt
#Date and time of measurement:
#Start: Tue Feb 16 2017, 23:14:06 UTC+1  
#Stop: Wed Feb 17 2017, 04:33:17 UTC+1

#Get Julian date of exposure
#JulianDay(month, day, year, hour)
t0 = JulianDay(2, 16, 2017, 22)
t1 = 14.0/(60*24) #Start time is 14 minutes past 22hrs
t2 = t1 + 5.31/24.0 #Total run time is 5.31 hours

Nvals = 10001
tvals = t0 + np.linspace(t1, t2, Nvals)

gammavals = np.zeros(Nvals)

#Calculate gamma from the LabVelocity
for i in range(Nvals):
    vs = -LabVelocity(tvals[i], lat_MPI, lon_MPI)
    vs_hat = vs/np.sqrt(np.sum(vs**2))
    rdet_hat = np.asarray([0,0,1.0])
    
    gammavals[i] = np.arccos(np.dot(vs_hat, rdet_hat))

#Load velocity distribution from file
def getVelDist(mstr, lsigstr, gamma_ind):
    Ngamvals = 11
    Nvvals = 61

    gamma_vals1, vvals1, fvals1 = np.loadtxt("../results/veldists/f_MPI_lmx" + mstr+"_lsig" + lsigstr + ".txt", unpack=True)
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
    integ = lambda x: eff(x)*((9.0/17.0)*dRdE(x, 27.0, m_x, sig, interpfun,vmax) + (8.0/17.0)*dRdE(x, 16.0, m_x, sig, interpfun,vmax))
    return quad(integ, E_min, E_max, epsrel=1e-4)[0]


#Calculate number of events as a function of gamma
Ne_list = np.zeros(11)
print "gamma    Ne"
for i in range(11):
    #We integrate from 1 -> 600 eV because the efficiency function takes care of the thresholds...
    Ne_list[i] = Nevents(1e-3, 600e-3, m_x, sigma_p, gamma_ind=i)*exposure

    print i*np.pi/10.0, Ne_list[i]

Ne_interp = interp1d(np.linspace(0, 1,11)*np.pi, Ne_list)

#Integrate over the values of gamma for a single day
#Note that the 5.31 hr exposure time is already included in the
#exposure, so need to correct for that here...
Ne_tot =  simps(Ne_interp(gammavals), tvals-t0)*(24.0/5.31)
print "Total number of events:", Ne_tot

#Append to number of events file
fname = "../results/Nevents/N_MPI_lmx" + '{0:.1f}'.format(np.log10(m_x)) + ".txt"
outarray = np.c_[np.log10(sigma_p), Ne_tot]

if (not os.path.isfile(fname)):
    htxt = "Total number of nucleus signal events at MPI (signal averaged over 1 day). Log10(m_x) = "+"{0:.1f}".format(np.log10(m_x))+"\nColumns: Log10(sigma/cm^2)      N_sig"
    np.savetxt(fname, outarray, header=htxt)
else:
    f_handle = file(fname, 'a')
    np.savetxt(f_handle, outarray)
    f_handle.close()


