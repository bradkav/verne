import sys
sys.path.append("../src/")

import numpy as np
import verne
from scipy.integrate import quad
from scipy.interpolate import interp1d

import matplotlib as mpl
font = { 'size'   : 16}
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


# Some random functions...
def q(A, E_R):
    mN = 1e3*A*0.9315 #in MeV
    
    return np.sqrt(2*mN*E_R*1e-3) #in MeV

def mu(mX, mA):
    return mX*mA*1.0/(mX+mA)

def ERmax(mX, mA, v):
    return 2*(mu(mX, mA)*v)**2/mA


m_x = 1

v_vals = np.linspace(0.1, 800, 100)

#Calculating correction factors...

print(">Calculating for Silicon...")
sig_corr_Si = v_vals*0.0
ER_corr_Si = v_vals*0.0
for i, v in enumerate(v_vals):
    sig_corr_Si[i] = quad(lambda x: verne.calcSIFormFactor(x*(1e6/(3e5*3e5))*ERmax(m_x, 0.9315*28, v), 28), 0, 1)[0]
    ER_corr_Si[i] = quad(lambda x: 2.0*x*verne.calcSIFormFactor(x*(1e6/(3e5*3e5))*ERmax(m_x, 0.9315*28, v), 28), 0, 1)[0]

print(">Calculating for Oxygen...")
sig_corr_Pb = v_vals*0.0
ER_corr_Pb = v_vals*0.0
for i, v in enumerate(v_vals):
    sig_corr_Pb[i] = quad(lambda x: verne.calcSIFormFactor(x*(1e6/(3e5*3e5))*ERmax(m_x, 0.9315*207, v), 207), 0, 1)[0]
    ER_corr_Pb[i] = quad(lambda x: 2.0*x*verne.calcSIFormFactor(x*(1e6/(3e5*3e5))*ERmax(m_x, 0.9315*207, v), 207), 0, 1)[0]

f, axarr = pl.subplots(2)


ax1 = axarr[0]
ax2 = axarr[1]

#Subplot 1 - cross section correction
ax1.plot(v_vals,sig_corr_Si, color='blue',label='Si', linewidth=1.5)
ax1.plot(v_vals,sig_corr_Pb, color='green', label='Pb', linewidth=1.5)
ax1.set_ylabel(r'Cross section correction, $\sigma(v)/\sigma(0)$', fontsize=12.0)
ax1.set_ylim(0, 1.1)
ax1.axhline(1.0, color='k', linestyle=':')

#Subplot 2 - recoil energy correction
ax2.plot(v_vals,ER_corr_Si, color='blue',label='Si', linewidth=1.5)
ax2.plot(v_vals,ER_corr_Pb, color='green', label='Pb', linewidth=1.5)

ax2.set_xlabel('v [km/s]', fontsize=20.0)
ax2.set_ylabel(r'$C_i(v) = <E_R>/E_\mathrm{max}$', fontsize=12.0)
ax2.set_ylim(0, 1.1)
ax2.axhline(1.0, color='k', linestyle=':')
ax2.legend(loc='best', frameon=False)

#f.subplots_adjust(hspace=0.1)
#pl.savefig("../plots/CorrectionFactors.pdf", bbox_inches="tight")

pl.show()
