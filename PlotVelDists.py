import numpy as np
import WIMpy.DMUtils as DMU
from numpy import pi
from scipy.integrate import quad, dblquad, trapz
import verne_analytic as VA
from verne_analytic import isoID
from skmonaco import mcquad
from timeit import default_timer as timer
from scipy.interpolate import interp2d, interp1d

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

print " Be careful - how much Iron is there actually..."

VA.loadIsotopes()

VA.loadPhiInterp()

rvals = np.linspace(0.0,650e4, 1000)


def getVelDist(lsigstr, vout, gamma):
    Nvals = 50
    Nthetavals = 10

    vvals, tvals, vfvals = np.loadtxt("results/VT/VT_lsig=" + lsigstr + ".txt", unpack=True)
    vigrid = vvals.reshape((Nthetavals,Nvals))
    tgrid = tvals.reshape((Nthetavals,Nvals))
    vfgrid = vfvals.reshape((Nthetavals,Nvals))

    est = 0.0*np.zeros((Nthetavals,len(vout)))
    for i in range(1,Nthetavals):
        theta = tgrid[i,0]
        #print theta
        #pl.figure()
        #print vfgrid[i,:]
        #pl.plot(vfgrid[i,:], vigrid[i,:])
        #pl.show()
        intp = interp1d(vfgrid[i,:], vigrid[i,:], bounds_error=False, fill_value=0.0)
        #pl.figure()
        #pl.plot(vout, intp(vout))
        #pl.show()
        
        dv = 1.0
        va = np.clip(vout+dv/2.0, 0, 800)
        vb = np.clip(vout-dv/2.0, 0, 800)
        
        deriv = np.abs(intp(va) - intp(vb))*1.0/(va-vb)
        vilist = intp(vout)
        
        est[i,:] = (deriv)*np.sin(theta)*(vilist**2)*VA.calcf_integ(vilist, theta, gamma)
        
    print tgrid[:,0]
    return trapz(est,tgrid[:,0], axis=0)



def getInterpFunc(lsigstr):
    Nvals = 50
    Nthetavals = 10

    vvals, tvals, vfvals = np.loadtxt("results/VT/VT_lsig=" + lsigstr + ".txt", unpack=True)
    vigrid = vvals.reshape((Nthetavals,Nvals))
    tgrid = tvals.reshape((Nthetavals,Nvals))
    vfgrid = vfvals.reshape((Nthetavals,Nvals))

    est = np.zeros(Nthetavals)
    for i in range(Nthetavals):
        intp = interp1d(vfgrid[i,:], vigrid[:,i])

    intp = interp2d(np.clip(vfgrid, 0.0, 800.0), tgrid, np.clip(vigrid, 0.0, 800.0))
    #print vigrid[-1,:],tgrid[-1,:], vfgrid[-1,:]
    vlist = np.linspace(0, 800.0, 100)
    tlist = np.linspace(np.pi/2, np.pi, 50)
    
    vlist2 = np.tile(vlist, 10)
    
    newgrid = np.zeros((100,50))
    for i,v in enumerate(vlist):
        for j,t in enumerate(tlist):
            newgrid[i,j] = intp(v,t)
    
    pl.figure()
    pl.contourf(vlist, tlist, newgrid, levels=np.linspace(0, 1000, 11))
    pl.colorbar()
    pl.show()

    return interp2d(np.clip(vfgrid, 0.0, 800.0), tgrid, np.clip(vigrid, 0.0, 800.0))

#print tgrid

v1 = np.linspace(0, 800, 100)

fv1 = getVelDist("-30.0", v1, np.pi*(140.0/180))
fv2 = getVelDist("-29.0", v1, np.pi*(140.0/180))
fv3 = getVelDist("-28.0", v1, np.pi*(140.0/180))

print fv1

#fv1 = v1*0.0
#fv2 = v1*0.0
#fv3 = v1*0.0

#int1 = getInterpFunc("-30.0")
#int2 = getInterpFunc("-29.0")
#int3 = getInterpFunc("-28.0")




#for i, v in enumerate(tqdm(v1)):
#    fv1[i] = VA.CalcFullF2(v, np.pi-0.1, int1)
#    fv2[i] = VA.CalcFullF2(v, np.pi-0.1, int2)
#    fv3[i] = VA.CalcFullF2(v, np.pi-0.1, int3)

#fv1 /= trapz(fv1, v1)
#fv2 /= trapz(fv2, v1)
#fv3 /= trapz(fv3, v1)

fig,ax1= pl.subplots(1)
ax1.plot(v1, DMU.calcf_SHM(v1),'k--',linewidth=1.5, label="Free")
ax1.plot(v1, fv1, linewidth=2.0, label=r'$\sigma_p = 10^{-30}$ cm$^2$')
ax1.plot(v1, fv2, linewidth=2.0, label=r'$\sigma_p = 10^{-29}$ cm$^2$')
ax1.plot(v1, fv3, linewidth=2.0, label=r'$\sigma_p = 10^{-28}$ cm$^2$')
ax1.axvline(DMU.vmin(10, 73, 1e5), linestyle=":", color='k')

ax1.set_xlabel(r'$v$ [km/s]')
ax1.set_ylabel(r'$f(v)$ [s/km]')
ax1.set_ylim(0, 5e-3)
pl.legend(loc='best')
pl.savefig('plots/SpeedDists.pdf', bbox_inches='tight')
pl.show()





