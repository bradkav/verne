import numpy as np
from numpy import pi
from scipy.integrate import quad, dblquad
import verne
from timeit import default_timer as timer
from matplotlib.ticker import MultipleLocator, FormatStrFormatter
import utils

from tqdm import tqdm

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


#15cm of Lead + 1cm of ancient lead (A = 207)
#Number density 3.3e22 per cc

#about 3cm copper A = 63.5
#Number density 8.5e22 per cc

#25cm of polyethylene
#Number density 

verne.loadIsotopes()

Nvals = 101

t1 = 1.0
theta = np.pi*t1
depth = 0.30
sigma_p = 2e-27
m_x = 1e5

verne.loadFFcorrections(m_x)

v_initial = np.linspace(0, 800, Nvals)
v_final_atmos = 0.0*v_initial
v_final_lab = 0.0*v_initial
v_final_full = 0.0*v_initial

for i, v in enumerate(v_initial):
    v_final_atmos[i] = verne.calcVfinal(v, theta,  depth, sigma_p, m_x, target="atmos")
    v_final_lab[i] = verne.calcVfinal(v_final_atmos[i], theta,  depth, sigma_p, m_x, target="earth")
    v_final_full[i] = verne.calcVfinal_shield_MPI(v_final_lab[i], sigma_p, m_x)
    

pl.figure()
pl.plot(v_initial, v_final_atmos,color='DarkBlue', linewidth=1.5, label='Atmos')
pl.plot(v_initial, v_final_lab,color='DarkOliveGreen', linewidth=1.5,label=' + Earth')
pl.plot(v_initial, v_final_full+2.0,color='DarkGoldenRod', linewidth=1.5, label=' + Shielding') #Add 2.0 so that it can be seen easily...
pl.plot([0, 800], [0, 800], 'k--')

#pl.title(r"$m_\chi = $" + utils.sciformat(m_x) + r"$, \sigma_p = $" + utils.sciformat(sigma_p) + r"$, \theta = $" + str(t1) + r" $\pi$")
pl.text(30,625, r"$m_\chi = $" + utils.sciformat(m_x) + r" $\mathrm{GeV}$" +\
         "\n" + r"$\sigma_p^{\mathrm{SI}} = $" + utils.sciformat(sigma_p) + r" $\mathrm{cm}^{2}$"+\
         "\nMPI (d = 0.3m)",\
         bbox=dict(boxstyle='round', facecolor='white', alpha=1.0) )

pl.text(475, 400, 'Atmos.', color='DarkBlue')
pl.text(455, 210, '+Earth', color='DarkOliveGreen')
pl.text(455, 75, '+Shield', color='DarkGoldenRod')

#pl.legend(loc='upper left', frameon=False)

pl.xlim(0, 800)
pl.ylim(0, 800)

pl.xlabel(r'$v_\mathrm{initial} \,[\mathrm{km/s}]$',fontsize=20.0)
pl.ylabel(r'$v_\mathrm{final} \,[\mathrm{km/s}]$', fontsize=20.0)

ax = pl.gca()
ax.xaxis.set_minor_locator(MultipleLocator(50))
ax.yaxis.set_minor_locator(MultipleLocator(50))

pl.savefig("plots/VelocityTransfer_nucleus.pdf", bbox_inches="tight")
pl.show()

