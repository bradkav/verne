import numpy as np
from LabFuncs import *

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

#SURF:
lat_SURF = 44.352 #N
lon_SURF = -103.751 #W

depth = 10.9e-3 #km
R_E = 6731.0 #km
r_det = R_E - depth


#NB: theta = 0 is from below, theta = 180 is from above
def calcEarthDepth(theta):
    trad = np.pi - theta*np.pi/180.0
    return -np.cos(trad)*r_det + np.sqrt((np.cos(trad)*r_det)**2 - (r_det**2 - R_E**2))
    
    
theta_vals = np.linspace(0, 180, 100)
    
pl.figure(figsize=(6,6))

pl.semilogy(theta_vals, calcEarthDepth(theta_vals), label='depth = 10.9m')

pl.xlabel(r"Incoming DM angle, $\gamma$")
pl.ylabel(r"Distance to detector [km]")

pl.legend(loc='best', frameon=False)
pl.xlim(0, 180)

#xlabs = [r'$'+str(x)+'^{\!\circ}$' for x in range(0, 200, 20)]
xlabs = [r'$0$', r'$\frac{\pi}{4}$',r'$\frac{\pi}{2}$',r'$\frac{3 \pi}{4}$',r'$\pi$']
a1 = pl.gca()
a1.set_xticks(np.linspace(0, 180, 5))
a1.set_xticklabels(xlabs)


pl.savefig("plots/depth_STAN.pdf", bbox_inches="tight")
pl.show()