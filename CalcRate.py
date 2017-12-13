import numpy as np
import WIMpy.DMUtils as DMU
from numpy import pi
from scipy.integrate import quad, dblquad, trapz
import verne
from LabFuncs import *
from matplotlib.ticker import MultipleLocator, FormatStrFormatter
import utils
from timeit import default_timer as timer
from scipy.interpolate import interp2d, interp1d
from matplotlib import cm
from tqdm import tqdm

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



def dfun(theta, phi):
    return np.arccos((np.sin(phi)*np.cos(theta) + np.cos(phi))/np.sqrt(2))*180/np.pi


def getGamma(t):
    vs = -LabVelocity(t, lat_STAN, lon_STAN)
    vs_hat = vs/np.sqrt(np.sum(vs**2))
    rdet_hat = np.asarray([0,0,1.0])
    
    return np.arccos(np.dot(vs_hat, rdet_hat))*180.0/np.pi


#SURF:
lat_SURF = 44.352 #N
lon_SURF = -103.751 #W

#Stanford
lat_STAN = +37.4 #N
lon_STAN = -122.2 #W

#NB: Soudan - where CDMS was IS A DIFFERENT PLACE!
lat_SOUDAN = 47.82

lat_LNGS = 42.454 #N
lon_LNGS = 13.576 #E

#Depth
l_STAN = 10.6e-3 #km
R_E = 6371.0 #km (MEAN RADIUS)

#(month, day, year, hour)
#t0 = JulianDay(1, 1, 2014, 1)
t0 = JulianDay(11, 1, 1998, 1)
Nvals = 10001


tvals = t0 + np.linspace(0, 365.0, Nvals)


t1 = 171.0
t2 = 172.0
tvals2 = t0 + np.linspace(t1, t2, Nvals)
#vs = np.zeros(3, Nvals)
gamma_list = np.zeros(Nvals)
alpha_list = np.zeros(Nvals)
beta_list = np.zeros(Nvals)

gamma_zoom = np.zeros(Nvals)

for i in tqdm(range(Nvals)):
    vs = -LabVelocity(tvals[i], lat_STAN, lon_STAN)
    vs_hat = vs/np.sqrt(np.sum(vs**2))
    
    vs2 = LabVelocity(tvals[i], 90.0, 0)
    
    rdet_hat = np.asarray([0,0,1.0])
    
    gamma_list[i] = np.arccos(np.dot(vs_hat, rdet_hat))*180.0/np.pi
    alpha_list[i] = np.arccos(vs2[2]/np.sqrt(np.sum(vs2**2)))*180.0/np.pi
    beta_list[i] = np.arctan2(vs2[0], vs2[1])*180.0/np.pi

for i in tqdm(range(Nvals)):
    vs = -LabVelocity(tvals2[i], lat_STAN, lon_STAN)
    vs_hat = vs/np.sqrt(np.sum(vs**2))
    
    vs2 = LabVelocity(tvals[i], 90.0, 0)
    
    rdet_hat = np.asarray([0,0,1.0])
    
    #gamma_zoom[i] = np.arccos(np.dot(vs_hat, rdet_hat))*180.0/np.pi
    gamma_zoom[i] = np.arccos(np.dot(vs_hat, rdet_hat))

dgamma_zoom = np.gradient(gamma_zoom, tvals2)

pl.figure()

pl.plot(tvals2, gamma_zoom/np.pi)
pl.show()

#print np.linspace(-29.1, -27.6, 16)/(1e0)

verne.loadIsotopes()

rvals = np.linspace(0.0,650e4, 1000)

v_th = DMU.vmin(20e-3, 27.0, 1e5)
print v_th

def getVelDist(lsigstr, gamma_ind):
    Ngamvals = 11
    Nvvals = 61
    
    rowvals = gamma_ind*61, 

    gamma_vals1, vvals1, fvals1 = np.loadtxt("results/veldists/f_lmx5.0_lsig" + lsigstr + ".txt", unpack=True)
    vvals = vvals1[gamma_ind*61:(gamma_ind+1)*61]
    fvals = fvals1[gamma_ind*61:(gamma_ind+1)*61]
    return vvals, fvals


def calcEta_final(v, interpfun):
    return quad(lambda x: interpfun(x)/x, v, 800.0)[0]

def dRdE(E, A, mx,sig,interpfun):  
    #sig = (1.973e-14*1.973e-14)*4.0*(DMU.reduced_m(1.0, mx))**2.0/np.pi

    int_factor = sig*DMU.calcSIFormFactor(E, A)*A**2
    
    return DMU.rate_prefactor(A, mx)*int_factor*calcEta_final(DMU.vmin(E, A, mx),  interpfun)

def Nevents(E_min, E_max, A, mx, sig, sigstr,gamma_ind = 10):
    vvals,fvals = getVelDist(sigstr, gamma_ind)
    interpfun = interp1d(vvals, fvals, kind='linear',bounds_error=False, fill_value=0.0)
    integ = lambda x: dRdE(x, A, mx, sig, interpfun)
    return quad(integ, E_min, E_max, epsrel=1e-4)[0]


Ne_list = np.zeros(11)

#print 10**(-30.6)

for i in range(11):
    #Ne_list[i] = Nevents(10, 100, 73, 1e5, 10**(-27.9), "-27.9", gamma_ind=i)*15.8
    Ne_list[i] = Nevents(20e-3, 10, 13, 1e5, 1e-30, "-30.0", gamma_ind=i)*(0.046e-3)
    print i, Ne_list[i]

Neinterp = interp1d(np.linspace(0, 1,11)*np.pi, Ne_list)

pl.figure()
pl.plot(tvals2-t0, Neinterp(gamma_zoom), label='Ne')
#pl.plot(tvals2, dgamma_zoom)
pl.legend()

pl.show()

print trapz(Neinterp(gamma_zoom), tvals2-t0)

#print Nevents(10, 100, 73, 1e5, 10**(-27.9), "-27.9", gamma_ind=6)*15.8 #kg days


