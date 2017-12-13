import numpy as np
import WIMpy.DMUtils as DMU
from numpy import pi
from scipy.integrate import quad, dblquad, trapz
import verne

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

#print np.linspace(-29.1, -27.6, 16)/(1e0)

verne.loadIsotopes()

rvals = np.linspace(0.0,650e4, 1000)

v_th = DMU.vmin(10.0, 73.0, m_x=1e5)

def getVelDist(lsigstr, gamma_ind):
    Ngamvals = 11
    Nvvals = 61
    
    rowvals = gamma_ind*61, 

    gamma_vals1, vvals1, fvals1 = np.loadtxt("results/veldists/f_lmx5.0_lsig" + lsigstr + ".txt", unpack=True)
    vvals = vvals1[gamma_ind*61:(gamma_ind+1)*61]
    fvals = fvals1[gamma_ind*61:(gamma_ind+1)*61]
    return vvals, fvals



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

v1 = np.linspace(0, 800, 100)

pl.figure()
ax1 = pl.gca()
slist = np.logspace(-30.6,-28.0, 14)
#print np.log10(slist)
slist = np.sort(np.append(slist, np.asarray([10**(-27.9), 10**(-27.8), 10**(-27.6)])))
strlist = ['{0:.1f}'.format(np.log10(sigma_p)) for sigma_p in slist]
#fname = "results/veldists/f_lmx" + '{0:.1f}'.format(np.log10(m_x)) + "_lsig" + '{0:.1f}'.format(np.log10(sigma_p)) + ".txt"

ax1.fill_between(np.linspace(0, v_th,100),0, 5,color='grey', alpha = 0.5, hatch="\\")

siglist = np.asarray([1e-30, 1e-29, 4e-29, 6e-29,1e-28, 1.25e-28])
cm_subsection = np.linspace(0.0, 0.85, len(siglist))
col_list = [ cm.Set1(x) for x in cm_subsection ]
ax1.plot(v1, DMU.calcf_SHM(v1),'k--',linewidth=1.5)
#for s in strlist:
for i,sig in enumerate(siglist):
    v, f = getVelDist("%.1f"%(np.log10(sig),), 7)
    labstr = " "
    #if (i == 0 or i == 10):
    #    labstr = str(i)
    ax1.plot(v, f, linewidth=2.0, color=col_list[i],label=str(int(sig*1e30)))
#ax1.plot(v1, fv2, linewidth=2.0, label=r'$\sigma_p = 10^{-29}$ cm$^2$')
#ax1.plot(v1, fv3, linewidth=2.0, label=r'$\sigma_p = 10^{-28}$ cm$^2$')
#ax1.axvline(DMU.vmin(10, 73, 1e5), linestyle=":", color='k')

ax1.set_xlabel(r'$v_f\, \,[\mathrm{km/s}]$',fontsize=20.0)
ax1.set_ylabel(r'$\tilde{f}(v_f) \,\,[\mathrm{s/km}]$',fontsize=20.0)
ax1.set_ylim(1e-7, 1e0)

ax1.yaxis.set_minor_locator(MultipleLocator(0.25))
ax1.xaxis.set_minor_locator(MultipleLocator(50))



m_x = 1e5
#sigma_p = 10**(-29.6)
#25,5*625/800.0
pl.text(30,4e-2, r"$m_\chi = $" + utils.sciformat(m_x) + r" $\mathrm{GeV}$" +\
         "\n" + r"$\gamma = 126^\circ$" + \
         "\nSUF (d = 10.6m)",\
         bbox=dict(boxstyle='round', facecolor='white', alpha=1.0) )

ax1.set_yscale("log")


#pl.text(575, 3e-1, r"$\sigma_p^\mathrm{SI} =$", fontsize=18.0)
pl.text(425, 3e-1, r"$\sigma_p^\mathrm{SI} = 10^{-30} \,\,\mathrm{cm}^2 \times$", fontsize=18.0)
#pl.text(610, 8e-4, r"$\times 10^{-30} \,\,\mathrm{cm}^2$", fontsize=18.0)
#pl.text(630, 4.55, r"above $(\gamma = \pi)$", fontsize=12.0)
#pl.text(630, 3.3, r"below $(\gamma = 0)$", fontsize=12.0)
pl.legend(loc='upper right',markerfirst=False,fontsize=14.0,frameon=False)
pl.savefig('plots/SpeedDists_xsec.pdf', bbox_inches='tight')
pl.show()





