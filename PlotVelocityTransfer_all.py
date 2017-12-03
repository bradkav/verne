import numpy as np
import WIMpy.DMUtils as DMU
from numpy import pi
from scipy.integrate import quad, dblquad
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


#15cm of Lead + 1cm of ancient lead (A = 207)
#Number density 3.3e22 per cc

#about 3cm copper A = 63.5
#Number density 8.5e22 per cc

#25cm of polyethylene
#Number density 

print "--------------"
print " CHECK THAT I'M NOT OUT BY A FACTOR OF 1e2 WITH LENGTHS!"
print " I'm sure I can speed this up too..."
print "--------------"

VA.loadIsotopes()
VA.loadPhiInterp()

Nvals = 10

t1 = 0.0
theta = np.pi*t1
depth = 10.6
sigma_p = 1e-30
m_x = 1e5


v_initial = np.linspace(0, 800, Nvals)
v_final = 0.0*v_initial
v_final_lab = 0.0*v_initial
v_final_full = 0.0*v_initial
calcVf = np.vectorize(VA.calcVfinal_full)


print VA.calcVfinal_full(600, theta,  10.6, sigma_p, m_x, target="earth")

for i, v in enumerate(v_initial):
    v_final[i] = VA.calcVfinal_full(v, theta,  depth, sigma_p, m_x, target="atmos")
    v_final_lab[i] = VA.calcVfinal_full(v, theta,  depth, sigma_p, m_x, target="full")
    v_final_full[i] = VA.calcVfinal_shield(v_final_lab[i], theta,  depth, sigma_p, m_x)
    

pl.figure()
pl.plot(v_initial, v_final,color='DarkBlue', linewidth=2.0, label='Atmos')
pl.plot(v_initial, v_final_lab,color='Brown', linewidth=2.0,label=' + Earth')
pl.plot(v_initial, v_final_full,color='DarkOrange', linewidth=2.0, label=' + Shielding')
pl.plot([0, 800], [0, 800], 'k--')

pl.title(r"$m_\chi = $" + utils.sciformat(m_x) + r"$, \sigma_p = $" + utils.sciformat(sigma_p) + r"$, \theta = $" + str(t1) + r" $\pi$")

#pl.text(600, 500, 'Atmosphere', color='b')

pl.legend(loc='upper left', frameon=False)

pl.xlim(0, 800)
pl.ylim(0, 800)

pl.xlabel(r'$v_\mathrm{initial}$ [km/s]',fontsize=18.0)
pl.ylabel(r'$v_\mathrm{final}$ [km/s]', fontsize=18.0)

pl.savefig("plots/VelocityTransfer.pdf", bbox_inches="tight")
pl.show()

