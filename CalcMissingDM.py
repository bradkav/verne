import numpy as np
import WIMpy.DMUtils as DMU
from numpy import pi
from scipy.integrate import quad, dblquad, quadrature, simps
from scipy.interpolate import interp1d
import verne

import MaxwellBoltzmann as MB

from timeit import default_timer as timer

import utils

from tqdm import tqdm

#Matplotlib ------------

import matplotlib as mpl
font = { 'size'   : 14}
mpl.rcParams['xtick.major.size'] = 5
mpl.rcParams['xtick.major.width'] = 1
mpl.rcParams['xtick.minor.size'] = 2
mpl.rcParams['xtick.minor.width'] = 1
mpl.rcParams['ytick.major.size'] = 5
mpl.rcParams['ytick.major.width'] = 1
mpl.rcParams['ytick.minor.size'] = 2
mpl.rcParams['ytick.minor.width'] = 1
mpl.rc('font', **font)

import matplotlib.pyplot as pl
#------------------------




verne.loadIsotopes()

Nvals = 1001

depth = 10.6
#sigma_p = 1.575e-28
sigma_p = 1e-29
m_x = 1e5

verne.loadFFcorrections(m_x)

gamma = np.pi
vesc = 533.0
v_e = np.sqrt(2)*156.0

tvals = np.linspace(0, 1, Nvals)
thetavals = np.pi*tvals


v_th = 1.0
#print verne.calcf_integ(220.0, 0, gamma), verne.calcf_integ(220.0, np.pi*0.5, gamma), verne.calcf_integ(220.0, np.pi, gamma)

#phi0 = 3*np.pi/4.0

#print verne.pathLength(0, np.pi), verne.pathLength(depth, np.pi)

#Varying phi doesn't do anything here because gamma = np.pi! --> Azimuthal symmetry


a = 1.0
b = 2*v_e*(-np.sin(gamma)*np.sin(np.pi-thetavals) + np.cos(gamma)*np.cos(np.pi-thetavals))
c = v_e**2 - vesc**2

v_initial = (-b + np.sqrt(b**2 - 4*a*c))/(2.0*a)
#v_initial2 = (-b - np.sqrt(b**2 - 4*a*c))/(2.0*a)

#pl.figure()
#pl.plot(tvals, v_initial1, label = "+")
#pl.plot(tvals, v_initial2, label = "-")
#pl.show()


v_final_lab = 0.0+0.0*v_initial
#print v_final_lab



#Now calculate the initial velocity which corresponds to v_cut as a final velocity



v_invert = 1.0*v_initial
for i in tqdm(range(Nvals-1, -1,-1)):
    v_invert[i] = verne.calcVinitial_full(v_th, thetavals[i],  depth, sigma_p, m_x, target="full")
    #print i, v_final_lab[i]
    if (v_invert[i] > v_initial[i]):
        v_invert[i] = v_initial[i]
        #break

pl.figure()
pl.plot(tvals, v_initial, 'k--')
pl.plot(tvals, v_invert,color='Brown', linewidth=2.0,label=' + Earth')
#pl.plot([0, 800], [0, 800], 'k--')

pl.title(r"$m_\chi = $" + utils.sciformat(m_x) + r"$, \sigma_p = $" + utils.sciformat(sigma_p) + r"$, \gamma = $" + str(gamma))

#pl.text(600, 500, 'Atmosphere', color='b')

#pl.legend(loc='upper left', frameon=False)

pl.xlim(0, 1.0)
pl.ylim(0, 800)

pl.xlabel(r'$\theta/\pi$',fontsize=18.0)
pl.ylabel(r'$v_\mathrm{initial}$ [km/s]', fontsize=18.0)


pl.show()

vmax = interp1d(thetavals, v_invert, kind='linear', bounds_error=False, fill_value=0.0)


#vcutoff = interp1d(thetavals, v_final_lab, kind='linear', bounds_error=False, fill_value=0)

tcinterp = interp1d(v_final_lab, thetavals, kind='linear', bounds_error=False, fill_value=0)


def finteg(v, theta):
    return v**2*np.sin(theta)*MB.calcf_integ(v, theta, gamma)

res = dblquad(finteg, 0, np.pi, lambda x: 0, lambda x: vmax(x), epsrel=1e-6)
print res

#vvals = np.linspace(0, 770, 100)
#fvals = 0.0*vvals
#for i,v in enumerate(vvals):
#    fvals[i] = v**2*quad(lambda x: np.sin(x)*verne.calcf_integ(v, x, gamma), 0,  tcinterp(vtest))[0]
#print "norm1 = ", simps(fvals, vvals)


