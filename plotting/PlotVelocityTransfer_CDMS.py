import sys
sys.path.append("../src/")

import numpy as np
import verne
from matplotlib.ticker import MultipleLocator, FormatStrFormatter
import utils

#Matplotlib ------------

import matplotlib as mpl
font = { 'size'   : 16, 'family': 'serif'}
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

verne.loadIsotopes()

Nvals = 101

t1 = 1.0
theta = np.pi*t1
depth = 10.6 #metres
sigma_p = 1e-28
m_x = 1e5

verne.loadFFcorrections(m_x)

v_initial = np.linspace(0, 800, Nvals)
v_final_atmos = 0.0*v_initial
v_final_lab = 0.0*v_initial
v_final_lab2 = 0.0*v_initial
v_final_lab3 = 0.0*v_initial
v_final_full = 0.0*v_initial

for i, v in enumerate(v_initial):
    v_final_atmos[i] = verne.calcVfinal(v, theta,  depth, sigma_p, m_x, target="atmos")
    v_final_lab[i] = verne.calcVfinal(v_final_atmos[i], theta,  depth, sigma_p, m_x, target="earth")
    v_final_lab2[i] = verne.calcVfinal(v_final_atmos[i], theta,  depth*1.04, sigma_p, m_x, target="earth")
    v_final_lab3[i] = verne.calcVfinal(v_final_atmos[i], theta,  depth*0.96, sigma_p, m_x, target="earth")
    v_final_full[i] = verne.calcVfinal_shield_SUF(v_final_lab[i], sigma_p, m_x)
    
    
print v_final_lab2/v_final_lab
print v_final_lab3/v_final_lab

pl.figure()
pl.plot(v_initial, v_final_atmos,color='DarkBlue', linewidth=1.5, label='Atmos')
pl.plot(v_initial, v_final_lab,color='DarkOliveGreen', linewidth=1.5,label=' + Earth')
pl.fill_between(v_initial, v_final_lab2,v_final_lab3,color='DarkOliveGreen', linewidth=1.0, alpha=0.5, zorder=0)
pl.plot(v_initial, v_final_full+2.0,color='DarkGoldenRod', linewidth=1.5, label=' + Shielding') #Add 2.0 so that it can be seen easily...
pl.plot([0, 800], [0, 800], 'k--')


pl.text(30,625, r"$m_\chi = $" + utils.sciformat(m_x) + r" $\mathrm{GeV}$" +\
         "\n" + r"$\sigma_p^{\mathrm{SI}} = $" + utils.sciformat(sigma_p) + r" $\mathrm{cm}^{2}$"+\
         "\nSUF (d = 10.6m)",\
         bbox=dict(boxstyle='round', facecolor='white', alpha=1.0) )

pl.text(475, 425, 'Atmos.', color='DarkBlue')
pl.text(455, 275, '+Earth', color='DarkOliveGreen')
pl.text(455, 75, '+Shield', color='DarkGoldenRod')


pl.xlim(0, 800)
pl.ylim(0, 800)

pl.xlabel(r'$v_\mathrm{initial} \,[\mathrm{km/s}]$',fontsize=20.0)
pl.ylabel(r'$v_\mathrm{final} \,[\mathrm{km/s}]$', fontsize=20.0)

ax = pl.gca()
ax.xaxis.set_minor_locator(MultipleLocator(50))
ax.yaxis.set_minor_locator(MultipleLocator(50))

pl.savefig("../plots/VelocityTransfer_CDMS.pdf", bbox_inches="tight")
pl.show()

