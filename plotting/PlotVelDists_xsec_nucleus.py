import sys
sys.path.append("../src/")


import numpy as np
from numpy import pi
import MaxwellBoltzmann as MB
from scipy.integrate import quad, trapz
from matplotlib.ticker import MultipleLocator, FormatStrFormatter
import utils
from scipy.interpolate import interp1d
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
v_th = MB.vmin(20e-3, 16.0, m_x=1e5)

#Load velocity distribution from file
def getVelDist(lsigstr, gamma_ind):
    Ngamvals = 11
    Nvvals = 61
    
    rowvals = gamma_ind*61, 

    gamma_vals1, vvals1, fvals1 = np.loadtxt("../results/veldists/f_MPI_lmx5.0_lsig" + lsigstr + ".txt", unpack=True)
    vvals = vvals1[gamma_ind*61:(gamma_ind+1)*61]
    fvals = fvals1[gamma_ind*61:(gamma_ind+1)*61]
    return vvals, fvals



v1 = np.linspace(0, 800, 1000)

pl.figure()
ax1 = pl.gca()

ax1.fill_between(np.linspace(0, v_th,100),0, 5,color='grey', alpha = 0.5, hatch="\\")

siglist = np.asarray([1e-28,6.3e-28, 1e-27,1.6e-27, 4e-27, 6e-27])
cm_subsection = np.linspace(0.0, 0.85, len(siglist))
col_list = [ cm.Set1(x) for x in cm_subsection ]
ax1.plot(v1, MB.calcf_SHM(v1),'k--',linewidth=1.5)

for i,sig in enumerate(siglist):
    v, f = getVelDist("%.2f"%(np.log10(sig),), 7)
    ax1.plot(v, f, linewidth=2.0, color=col_list[i],label=str(int(sig*1.0001e28)))

ax1.set_xlabel(r'$v_f\, \,[\mathrm{km/s}]$',fontsize=20.0)
ax1.set_ylabel(r'$\tilde{f}(v_f) \,\,[\mathrm{s/km}]$',fontsize=20.0)
ax1.set_ylim(1e-7, 1e0)

ax1.yaxis.set_minor_locator(MultipleLocator(0.25))
ax1.xaxis.set_minor_locator(MultipleLocator(50))

pl.text(30,4e-2, r"$m_\chi = $" + utils.sciformat(1e5) + r" $\mathrm{GeV}$" +\
         "\n" + r"$\gamma = 126^\circ$" + \
         "\nMPI (d = 0.3m)",\
         bbox=dict(boxstyle='round', facecolor='white', alpha=1.0) )

ax1.set_yscale("log")

pl.text(425, 3e-1, r"$\sigma_p^\mathrm{SI} = 10^{-28} \,\,\mathrm{cm}^2 \times$", fontsize=18.0)
pl.legend(loc='upper right',markerfirst=False,fontsize=14.0,frameon=False)
pl.savefig('../plots/SpeedDists_xsec_nucleus.pdf', bbox_inches='tight')
pl.show()





