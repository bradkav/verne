import sys
sys.path.append("../src/")

import numpy as np
import MaxwellBoltzmann as MB
from numpy import pi
from scipy.integrate import quad, trapz
from matplotlib.ticker import MultipleLocator, FormatStrFormatter
import utils
from scipy.interpolate import interp2d, interp1d
from matplotlib import cm

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


#Threshold velocity
v_th = MB.vmin(10.0, 73.0, m_x=1e5)

#Load velocity distribution from file
def getVelDist(lsigstr, gamma_ind):
    Ngamvals = 11
    Nvvals = 61
    
    rowvals = gamma_ind*61, 

    #gamma_vals1, vvals1, fvals1 = np.loadtxt("../results/veldists/f_SUF_lmx5.0_lsig" + lsigstr + ".txt", unpack=True)
    gamma_vals1, vvals1, fvals1 = np.loadtxt("../results/veldists/f_SI_full_mx100000.000_lsig" + lsigstr + ".txt", unpack=True)
    vvals = vvals1[gamma_ind*61:(gamma_ind+1)*61]
    fvals = fvals1[gamma_ind*61:(gamma_ind+1)*61]
    return vvals, fvals



v1 = np.linspace(0, 800, 100)

pl.figure()
ax1 = pl.gca()

cm_subsection = np.linspace(0.0, 1.0, 11)
col_list = [ cm.viridis(x) for x in cm_subsection ]
ax1.plot(v1, 1e3*MB.calcf_SHM(v1),'k--',linewidth=1.5)
s = "-29.60"

ax1.fill_between(np.linspace(0, v_th,100),0, 5,color='grey', alpha = 0.5, hatch="\\")

for i in range(10,-1, -1):
    v, f = getVelDist(s, i)
    ax1.plot(v, 1e3*f, linewidth=2.0, color=col_list[i], label=" ")

ax1.set_xlabel(r'$v_f\, \,[\mathrm{km/s}]$',fontsize=20.0)
ax1.set_ylabel(r'$\tilde{f}(v_f) \,\,[10^{-3} \,\mathrm{s/km}]$',fontsize=20.0)
ax1.set_ylim(0, 5)

ax1.yaxis.set_minor_locator(MultipleLocator(0.25))
ax1.xaxis.set_minor_locator(MultipleLocator(50))



m_x = 1e5
sigma_p = 10**(-29.6)

pl.text(30,5*625/800.0, r"$m_\chi = $" + utils.sciformat(m_x) + r" $\mathrm{GeV}$" +\
         "\n" + r"$\sigma_p^{\mathrm{SI}} = $" + utils.sciformat_1(sigma_p) + r" $\mathrm{cm}^{2}$" + \
         "\nSUF (d = 10.6m)",\
         bbox=dict(boxstyle='round', facecolor='white', alpha=1.0) )


pl.text(375, 3.9, r"Average DM flux from...", fontsize=12.0)
pl.text(610, 4.55, r"above $(\gamma = 180^\circ)$", fontsize=12.0)
pl.text(620, 3.3, r"below $(\gamma = 0^\circ)$", fontsize=12.0)
pl.legend(loc=[0.8, 0.70], fontsize=6.0,labelspacing=0.001, handlelength=10.0, frameon=False)
pl.savefig('../plots/SpeedDists_gamma.pdf', bbox_inches='tight')
pl.show()





