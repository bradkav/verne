import numpy as np
import WIMpy.DMUtils as DMU
from numpy import pi
from scipy.integrate import quad, dblquad
import verne_analytic as VA
from verne_analytic import isoID
from skmonaco import mcquad
from timeit import default_timer as timer

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

print " Be careful - how much Iron is there actually..."

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

VA.LoadVelocityFunctions()

Fv, Fv_inv = VA.calcVelocityFunction(isoID["Fe"])

v_vals = np.linspace(1e-2, 800, 1000)

start = timer()
result, error = mcquad(lambda x: x[0]**2*np.sin(x[1])*VA.calcf_integ(x[0], x[1],0.0), npoints=1e5, xl=[0., 0.], xu=[760.0, np.pi])
end = timer()
print " Time taken: ", end - start
print "{} +/- {}".format(result,error)

#pl.figure()
#pl.plot(v_vals, Fv(v_vals))
#pl.show()

VA.calcDeff_interp(10.6, isoID["Fe"])

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
Nvals = 10

vlist = np.linspace(0, 750, Nvals+1)
etalist = np.zeros(Nvals)
errorlist = np.zeros(Nvals)
for i in range(Nvals):
    start = timer()
    etalist[i],errorlist[i] = VA.CalcFullEta(vlist[i], np.pi, 1e-28, 1e5, isoID["Fe"], vb=vlist[i+1])
    end = timer()
    print " ", i+1, "- Time taken: ", end - start

print etalist
print errorlist

fulleta = np.sum(etalist) - np.cumsum(etalist)
print fulleta

fig,ax1= pl.subplots(1)
ax1.plot(vlist[1:], fulleta, 'g^')
#ax1.errorbar(vlist[1:], fulleta, yerr=errorlist, fmt='^', color='g')
ax1.plot(np.linspace(0,750,100), DMU.calcEta(np.linspace(0,750,100)), 'b-')

ax1.set_yscale('log')


ax1.set_xlabel(r'$v_\mathrm{min}$ [km/s]')
ax1.set_ylabel(r'$\eta(v_\mathrm{min})$ [s/km]')

ax1.set_xlim(0, 750)
ax1.set_ylim(1e-10, 1e0)

pl.show()

