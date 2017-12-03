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

#NB: Soudan - where CDMS was IS A DIFFERENT PLACE!
lat_SOUDAN = 47.82

lat_LNGS = 42.454 #N
lon_LNGS = 13.576 #E

t0 = JulianDay(1, 1, 2014, 1)

Nvals = 366
tvals = t0 + np.linspace(0, 365, Nvals)
#vs = np.zeros(3, Nvals)
alpha_list = np.zeros(Nvals) 

for i in range(Nvals):
    vs = LabVelocity(tvals[i], 90.0, 0)
    alpha_list[i] = np.arccos(vs[2]/np.sqrt(np.sum(vs**2)))*180.0/np.pi
    
    
    
pl.figure()
pl.xlabel("Days from 1 Jan 2014")
pl.ylabel(r"Latitude of mean incoming DM, $\alpha$")
pl.plot(tvals-t0, alpha_list)
#The alpha and lambda angles are complementary!
pl.axhline(90-lat_SURF, linestyle="--", color="k")
pl.text(300, 90-lat_SURF+0.25, "SURF")

pl.xlim(0, 365)
pl.savefig("plots/Alpha.pdf", bbox_inches="tight")
pl.show()