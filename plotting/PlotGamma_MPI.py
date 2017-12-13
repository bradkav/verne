import sys
sys.path.append("../src/")

import numpy as np
from LabFuncs import *

import matplotlib as mpl
import matplotlib.patches as patches
font = { 'size'   : 16, 'family':'serif'}
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


def dfun(theta, phi):
    return np.arccos((np.sin(phi)*np.cos(theta) + np.cos(phi))/np.sqrt(2))*180/np.pi


def getGamma(t):
    vs = -LabVelocity(t, lat_MPI, lon_MPI)
    vs_hat = vs/np.sqrt(np.sum(vs**2))
    rdet_hat = np.asarray([0,0,1.0])
    
    return np.arccos(np.dot(vs_hat, rdet_hat))*180.0/np.pi



#MPI Munich
lat_MPI = +48.1 #N
lon_MPI = +11.57 #W

#(month, day, year, hour)
#t0 = JulianDay(1, 1, 2014, 1)
t0 = JulianDay(1, 1, 2017, 1)
Nvals = 10001


tvals = t0 + np.linspace(0, 365.0, Nvals)


t1 = 170.0
t2 = 190.0
tvals2 = t0 + np.linspace(t1, t2, Nvals)
#vs = np.zeros(3, Nvals)
gamma_list = np.zeros(Nvals)
alpha_list = np.zeros(Nvals)
beta_list = np.zeros(Nvals)

gamma_zoom = np.zeros(Nvals)

for i in range(Nvals):
    vs = -LabVelocity(tvals[i], lat_MPI, lon_MPI)
    vs_hat = vs/np.sqrt(np.sum(vs**2))
    
    vs2 = LabVelocity(tvals[i], 90.0, 0)
    
    rdet_hat = np.asarray([0,0,1.0])
    
    gamma_list[i] = np.arccos(np.dot(vs_hat, rdet_hat))*180.0/np.pi
    alpha_list[i] = np.arccos(vs2[2]/np.sqrt(np.sum(vs2**2)))*180.0/np.pi
    beta_list[i] = np.arctan2(vs2[0], vs2[1])*180.0/np.pi

for i in range(Nvals):
    vs = -LabVelocity(tvals2[i], lat_MPI, lon_MPI)
    vs_hat = vs/np.sqrt(np.sum(vs**2))
    
    vs2 = LabVelocity(tvals[i], 90.0, 0)
    
    rdet_hat = np.asarray([0,0,1.0])
    
    gamma_zoom[i] = np.arccos(np.dot(vs_hat, rdet_hat))*180.0/np.pi



print " Maximum gamma: ", np.max(gamma_list)
print " Time of max: ", tvals[np.argmax(gamma_list)]-t0 
maxind = np.argmax(gamma_list)
    
    
pl.figure()
pl.xlabel("Days from 1 Jan 2017")
pl.ylabel(r"$\gamma = \cos^{-1} (\langle \hat{\mathbf{v}}_\chi \rangle\cdot \hat{\mathbf{r}}_\mathrm{det})$")

pl.fill_between(tvals-t0,(180-alpha_list)-(90-lat_MPI),180.0-np.abs(alpha_list - (90-lat_MPI)), color='b', alpha=0.4)
pl.plot(tvals-t0, (180-alpha_list)-(90-lat_MPI), 'b-')
pl.plot(tvals-t0, 180.0-np.abs(alpha_list - (90-lat_MPI)), 'b-')

pl.axhline(180-(90-lat_MPI), linestyle="--", color="k")


pl.xlim(0,365)
pl.ylim(0,180)

a1 = pl.gca()
a1.add_patch(
     patches.Rectangle(
        (170.0, 70.0),
        20.0,
        110.0,
        fill=False,      # remove background
        linewidth=0.75
     ))
ylabs = [r'$'+str(x)+'^{\!\circ}$' for x in range(0, 200, 20)]
a1.set_yticklabels(ylabs)
a1.text(5, 172.6, "DM mostly from above",fontsize=10.0)
a1.text(5, 5,"DM mostly from below", fontsize=10.0)


pl.plot([190.0, 336.5], [180.0, 81.0], 'k-', linewidth=0.75)
pl.plot([170.0, 195.0], [70.0, 14.0], 'k-', linewidth=0.75)

a = pl.axes([.54, .16, .3, .3])
a.plot(tvals2-t0, gamma_zoom)
a.yaxis.tick_right()
pl.xlim(170, 190)
pl.ylim(70, 180)

a.set_xticklabels([170, 175, 180, 185, 190], fontsize='small')

a.set_yticks([80, 100, 120, 140, 160, 180])
a.set_yticklabels([r'$'+str(x)+'^{\!\circ}$' for x in range(80, 200, 20)], fontsize='small')

pl.savefig("../plots/Gamma_MPI.pdf",bbox_inches="tight")
pl.show()