import numpy as np
import WIMpy.DMUtils as DMU
from numpy import pi
from scipy.integrate import quad, dblquad, trapz
import verne_analytic as VA
from verne_analytic import isoID
from skmonaco import mcquad
from timeit import default_timer as timer

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

VA.loadIsotopes()
VA.loadPhiInterp()

rvals = np.linspace(0.0,650e4, 1000)

#pl.figure()
#pl.plot(rvals/637.1e4, VA.dens_interp[isoID["Fe"]](rvals))
#pl.show()

tvals = np.linspace(0, np.pi,100)
#pl.figure()
#pl.semilogy(tvals, VA.pathLength(1500, tvals))
#pl.show()

Deff = tvals*0.0
Deff2 = tvals*0.0

#for i, t in enumerate(tvals):
#    Deff[i] = VA.calcDeff(1500, t, isoID["Fe"])
    #Deff2[i] = VA.calcDeff(1500, t, isoID["O"])

#pl.figure()
#pl.semilogy(tvals/np.pi, Deff)
#pl.semilogy(tvals/np.pi, Deff2)
#pl.show()

#print VA.effectiveXS(1e-28, 1e5, isoID["Fe"])

Vinit = np.vectorize(VA.calcVinitial)

eta_int = np.vectorize(VA.eta_integrand)

#pl.figure()
#pl.plot(tvals,Vinit(220, tvals, 1e-28, 1e5,  isoID["Fe"]))
#pl.show()

#vg, tg = np.meshgrid(v_vals, tvals, indexing='ij')
#v_i = Vinit(vg, tg, 1e-29, 1e5,  isoID["Fe"])
#etagrid = eta_int(vg, tg, 0, 1e-30, 1e5,  isoID["Fe"])

#pl.figure()
#pl.contourf(vg,tg,etagrid)
#pl.colorbar()
#pl.show()


"""
pl.figure()
pl.plot(v_vals,VA.calcVinitial(v_vals, np.pi/2.0+0.1, 1e-28, 1e5,  isoID["Fe"]))
pl.show()

pl.figure()
pl.plot(v_vals,VA.calcVfinal(v_vals, np.pi, 1e-28, 1e5,  isoID["Fe"]))
pl.show()
"""
print " Be careful about the ranges of my interpolation functions..."
print " Must be a problem with my interpolation functions..."

#start = timer()
#print VA.eta_integrand(200, 0.1,  0, 1e-30, 1e5, isoID["Fe"])
#end = timer()
#print " Time taken: ", end - start
mx = 1e5

#siglist = [-29.0, -28.0, -27.0]
siglist = [-31.0, -30.0, -29.0]
sigvals = 10**np.asarray(siglist)
labs = [r'$\sigma_p = 10^{'+str(int(lsig))+'}$ cm$^2$' for lsig in siglist]
#labs = [r'$s = 1$',r'$s = 10$', r'$s = 100$',r'$s = 1000$']

Nvals = 5


vlist = np.logspace(np.log10(1), np.log10(750), Nvals)
print vlist
flist = np.zeros((3, Nvals))
for j in range(3):
    for i in tqdm(range(Nvals)):
        flist[j,i] = VA.CalcFullF_full(vlist[i], np.pi, 10.6, sigvals[j], mx, "earth", np.pi/2.0)
        #print vlist[i], flist[j,i]
    #print " I'm renormalising by fudging it!!!"
    #norm1 = -np.trapz(vlist, flist[j,:])
    #print norm1
    #flist[j,:] /= norm1

#print flist[2,:]

fig,ax1= pl.subplots(1, figsize=(8,6))
#ax1.errorbar(vlist[1:], fulleta, yerr=errorlist, fmt='^', color='g')
ax1.plot(np.linspace(0,750,100), DMU.calcf_SHM(np.linspace(0,750,100)),'k--',linewidth=1.0, label="Free")
for j in range(3):
    print trapz(flist[j,1:],vlist[1:])
    ax1.plot(vlist, flist[j,:],label=labs[j], linewidth=1.5)
ax1.set_yscale('log')

ax1.axvline(DMU.vmin(10, 73, 1e5), linestyle=":", color='k')

ax1.set_title(r'$m_\chi = 10^5$ GeV')

ax1.set_xlabel(r'$v$ [km/s]')
ax1.set_ylabel(r'$f(v)$ [s/km]')

#ax1.set_xlim(0, 750)
ax1.set_ylim(1e-8, 5e-3)

pl.legend(loc='best', fontsize=14.0)
pl.savefig('plots/sigma_STAN_log.pdf', bbox_inches="tight")
pl.show()

