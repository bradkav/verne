import numpy as np
import WIMpy.DMUtils as DMU
from numpy import pi
from scipy.integrate import quad, dblquad, quadrature, simps
from scipy.interpolate import interp1d
import verne_analytic as VA
from verne_analytic import isoID
from skmonaco import mcquad
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




VA.loadIsotopes()
VA.loadPhiInterp()

Nvals = 1000

depth = 10.6
sigma_p = 1e-28
m_x = 1e5

gamma = np.pi
vesc = 533.0
v_e = np.sqrt(2)*156.0

tvals = np.linspace(0, 1, Nvals)
thetavals = np.pi*tvals


v_th = 1.0
#print VA.calcf_integ(220.0, 0, gamma), VA.calcf_integ(220.0, np.pi*0.5, gamma), VA.calcf_integ(220.0, np.pi, gamma)

#phi0 = 3*np.pi/4.0

#print VA.pathLength(0, np.pi), VA.pathLength(depth, np.pi)

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

for i in tqdm(range(Nvals-1, -1,-1)):
    
    v_final_lab[i] = VA.calcVfinal_full(v_initial[i], thetavals[i],  depth, sigma_p, m_x, target="earth")
    #print i, v_final_lab[i]
    if (v_final_lab[i] < v_th):
        break
    
tcut = np.min(thetavals[v_final_lab > v_th]) + thetavals[1] - thetavals[0]
vcut = np.max(v_final_lab)
print tcut
print vcut

pl.figure()
pl.plot(tvals, v_initial, 'k--')
pl.plot(tvals, v_final_lab,color='Brown', linewidth=2.0,label=' + Earth')
#pl.plot([0, 800], [0, 800], 'k--')

pl.title(r"$m_\chi = $" + utils.sciformat(m_x) + r"$, \sigma_p = $" + utils.sciformat(sigma_p) + r"$, \gamma = $" + str(gamma))

#pl.text(600, 500, 'Atmosphere', color='b')

#pl.legend(loc='upper left', frameon=False)

pl.xlim(0, 1.0)
pl.ylim(0, 800)

pl.xlabel(r'$\theta/\pi$',fontsize=18.0)
pl.ylabel(r'$v_\mathrm{final}$ [km/s]', fontsize=18.0)

pl.savefig("plots/VelocityCutOff.pdf", bbox_inches="tight")
pl.show()

#vcutoff = interp1d(thetavals, v_final_lab, kind='linear', bounds_error=False, fill_value=0)

tcinterp = interp1d(v_final_lab, thetavals, kind='linear', bounds_error=False, fill_value=0)

"""
v1 = np.linspace(0, 800, 100)
pl.figure()
pl.plot(v_final_lab, thetavals, 'b+')
pl.plot(v1, tcinterp(v1))
pl.show()
"""

start = timer()

vtest = 85
#vtest = 100
delta = np.pi - tcinterp(vtest)
tlist = tcinterp(vtest)+np.append(0.0, np.logspace(-4, 0, 100)*delta)
fint = tlist*0.0
for i in tqdm(range(len(tlist))):
    fint[i] = VA.f_integrand_full(vtest, tlist[i], gamma, depth, sigma_p, m_x, target="earth")

end = timer()
t1 = end - start

print simps(fint, tlist)
pl.figure()
pl.plot(tlist/np.pi, fint, '+')
pl.axvline(tcinterp(vtest)/np.pi, linestyle='--', color='k')
pl.show()

print "WHICH ANSWER IS CORRECT!"
print "DO SOME REFINEMENT!"



print simps(fint, tlist), " (time taken:", t1, ")"

start = timer()
other = VA.CalcFullF_full(vtest, gamma, depth, sigma_p, m_x, "earth", tcut = tcinterp(vtest))
end = timer()
t2 = end - start

print other, " (time taken:", t2, ")"

pl.figure()
pl.plot(tlist/np.pi, fint)
pl.axvline(tcinterp(vtest)/np.pi, linestyle='--', color='k')
pl.show()



vvals = np.linspace(0, 770, 50)
fvals = 0.0*vvals
for i,v in enumerate(vvals):
    fvals[i] = v**2*quad(lambda x: np.sin(x)*VA.calcf_integ(v, x, gamma), 0,  tcinterp(vtest))[0]
print "norm1 = ", np.trapz(fvals, vvals)



vlist = np.logspace(np.log10(v_th), np.log10(vcut), 50)
fout = 0.0*vlist
fnorm = 0.0*vlist
for i in tqdm(range(len(vlist))):
    #if (vlist[i] > vcut):
    #    fout[i] =  0
    #else:
    #print tcinterp(vlist[i])/np.pi
    #print tcut/np.pi
    print "Calculating for v =", vlist[i]
    fnorm[i] = vlist[i]**2*quad(lambda t: np.sin(t)*VA.calcf_integ(vlist[i], t, gamma), 0, np.pi)[0]
    fout[i] = VA.CalcFullF_full(vlist[i], gamma, depth, sigma_p, m_x, "earth", tcut = tcinterp(vlist[i]))

#integ1 = lambda x: VA.CalcFullF_full(x[0], gamma, depth, sigma_p, m_x, "earth", tcut = tcinterp(x[0]))

#Norm = quad(integ1, 0, vcut, rtol=1e-2, maxiter=20)[0] 
#print Norm
print np.trapz(fout, vlist)

pl.figure()
pl.plot(vlist, fout)
pl.plot(vlist, fnorm)

pl.figure()
pl.semilogy(vlist, fout)
pl.semilogy(vlist, fnorm)
pl.show()

