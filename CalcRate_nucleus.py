import numpy as np
import WIMpy.DMUtils as DMU
from numpy import pi
from scipy.integrate import quad, dblquad, trapz, simps
import verne
from LabFuncs import *
from matplotlib.ticker import MultipleLocator, FormatStrFormatter
import utils
from scipy.special import erf
from timeit import default_timer as timer
from scipy.interpolate import interp2d, interp1d
from matplotlib import cm
from tqdm import tqdm
import argparse
import os.path


#Matplotlib ------------

import matplotlib as mpl
font = { 'size'   : 16, 'family':'serif'}
mpl.rcParams['xtick.major.size'] = 5
mpl.rcParams['xtick.major.width'] = 1
mpl.rcParams['xtick.minor.size'] = 3
mpl.rcParams['xtick.minor.width'] = 1
mpl.rcParams['ytick.major.size'] = 5
mpl.rcParams['ytick.major.width'] = 1
mpl.rcParams['ytick.minor.size'] = 3
mpl.rcParams['ytick.minor.width'] = 1
mpl.rc('font', **font)

import matplotlib.pyplot as pl
#------------------------

#Parse the arguments!
parser = argparse.ArgumentParser(description='...')
parser.add_argument('-m_x','--m_x', help='DM mass in GeV', type=float,default = 1e5)
parser.add_argument('-sigma_p','--sigma_p', help='DM-nucleon cross section, sigma_p in cm^2', type=float, required=True)
args = parser.parse_args()
m_x = args.m_x
sigma_p = args.sigma_p



def eff(E):
    E_th = 19.7e-3
    sig_th = 3.83e-3
    return 0.5*(1+erf((E - E_th)/(np.sqrt(2)*sig_th)))


#MPI Munich
lat_MPI = +48.1 #N
lon_MPI = +11.57 #W

#(month, day, year, hour)
#t0 = JulianDay(1, 1, 2014, 1)
t0 = JulianDay(1, 1, 2017, 1)
Nvals = 10001

#June 20th...
t1 = 171.0
t2 = 172.0
tvals = t0 + np.linspace(t1, t2, Nvals)

gammavals = np.zeros(Nvals)

for i in range(Nvals):
    vs = -LabVelocity(tvals[i], lat_MPI, lon_MPI)
    vs_hat = vs/np.sqrt(np.sum(vs**2))
    rdet_hat = np.asarray([0,0,1.0])
    
    #gamma_zoom[i] = np.arccos(np.dot(vs_hat, rdet_hat))*180.0/np.pi
    gammavals[i] = np.arccos(np.dot(vs_hat, rdet_hat))


def getVelDist(mstr, lsigstr, gamma_ind):
    Ngamvals = 11
    Nvvals = 61

    gamma_vals1, vvals1, fvals1 = np.loadtxt("results/veldists/f_MPI_lmx" + mstr+"_lsig" + lsigstr + ".txt", unpack=True)
    vvals = vvals1[gamma_ind*61:(gamma_ind+1)*61]
    fvals = fvals1[gamma_ind*61:(gamma_ind+1)*61]
    return vvals, fvals


def calcEta_final(v, interpfun, vmax):
    return quad(lambda x: interpfun(x)/x, v, vmax*1.1)[0]

def dRdE(E, A, mx,sig,interpfun,vmax):  
    #sig = (1.973e-14*1.973e-14)*4.0*(DMU.reduced_m(1.0, mx))**2.0/np.pi

    int_factor = sig*DMU.calcSIFormFactor(E, A)*A**2
    
    return DMU.rate_prefactor(A, mx)*int_factor*calcEta_final(DMU.vmin(E, A, mx),  interpfun, vmax)

def Nevents(E_min, E_max, m_x, sig, gamma_ind = 10):
    sigstring = '{0:.2f}'.format(np.log10(sig))
    mstring = '{0:.1f}'.format(np.log10(m_x))
    
    vvals,fvals = getVelDist(mstring, sigstring, gamma_ind)
    vmax = np.max(vvals)
    #pl.figure()
    #pl.plot(vvals, fvals)
    #pl.show()
    interpfun = interp1d(vvals, fvals, kind='linear',bounds_error=False, fill_value=0.0)
    integ = lambda x: eff(x)*((9.0/17.0)*dRdE(x, 27.0, m_x, sig, interpfun,vmax) + (8.0/17.0)*dRdE(x, 16.0, m_x, sig, interpfun,vmax))
    return quad(integ, E_min, E_max, epsrel=1e-4)[0]


exposure = 0.046e-3

Ne_list = np.zeros(11)

for i in range(11):
    #Ne_list[i] = Nevents(10, 100, 73, 1e5, 10**(-27.9), "-27.9", gamma_ind=i)*15.8
    Ne_list[i] = Nevents(1e-3, 600e-3, m_x, sigma_p, gamma_ind=i)*exposure
    print i, Ne_list[i]

Ne_interp = interp1d(np.linspace(0, 1,11)*np.pi, Ne_list)

Ne_tot =  simps(Ne_interp(gammavals), tvals-t0)
print "Total number of events:", Ne_tot

fname = "results/Nevents/N_MPI_lmx" + '{0:.1f}'.format(np.log10(m_x)) + ".txt"

outarray = np.c_[np.log10(sigma_p), Ne_tot]

if (not os.path.isfile(fname)):
    htxt = "Total number of nucleus signal events at MPI (signal averaged over 1 day). Log10(m_x) = "+"{0:.1f}".format(np.log10(m_x))+"\nColumns: Log10(sigma/cm^2)      N_sig"
    np.savetxt(fname, outarray, header=htxt)
else:
    f_handle = file(fname, 'a')
    #f_handle.write(b'\n')
    #np.savetxt(f_handle, "\n")
    np.savetxt(f_handle, outarray)
    f_handle.close()

#print Nevents(10, 100, 73, 1e5, 10**(-27.9), "-27.9", gamma_ind=6)*15.8 #kg days


