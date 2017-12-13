import numpy as np
from LabFuncs import *
import WIMpy.DMUtils as DMU
from scipy.integrate import quad
from scipy.interpolate import interp1d

from tqdm import tqdm

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

print 1e6*ERmax(1e10, 0.9315*16, 750/3.0e5)
print 1e6*ERmax(1e10, 0.9315*56, 750/3.0e5)

m_x = 1e5

v_vals = np.linspace(0.1, 800, 100)

#Calculating correction factors...

#v_vals_log = np.logspace(np.log10(1e-2), np.log10(800), 200)
A0 = 56
corr_fact = v_vals*0.0
for i, v in enumerate(v_vals):
        corr_fact[i] = quad(lambda x: 2.0*x*DMU.calcSIFormFactor(x*(1e6/(3e5*3e5))*ERmax(m_x, 0.9315*A0, v), A0), 0, 1)[0]
corr_interp = interp1d(v_vals, corr_fact, kind='linear')

#Calculating velocity factors...

velfunc = v_vals*0.0
for i in tqdm(range(len(v_vals))):
    velfunc[i] = quad(lambda x: 1.0/(x*corr_interp(x)), 0.1, v_vals[i])[0]
    #print velfunc[i]
velinterp = interp1d(v_vals, velfunc, kind='linear')
velinterp_inv = interp1d(velfunc, v_vals, kind='linear', bounds_error=False, fill_value = 1e-10)

print -1.0/1.0 + velinterp(750)
print velinterp_inv(-1.0/1.0 + velinterp(750))

pl.figure()
pl.plot(v_vals, velfunc)
pl.show()



def finalV(v0, d, lam=1.0, kind="full"):
    if (kind == "full"):
        return velinterp_inv(-d/lam + velinterp(v0))
    if (kind == "approx"):
        return v0*np.exp(-d/lam)
    
pl.figure()
pl.plot(v_vals, finalV(v_vals, d=0.5, lam=1.0), 'b-',label="d = 0.5")
pl.plot(v_vals, finalV(v_vals, d=2.0, lam=1.0), 'g-',label="d = 2")
pl.plot(v_vals, finalV(v_vals, d=5.0, lam=1.0), 'r-',label="d = 5")
pl.plot(v_vals, finalV(v_vals, d=20.0, lam=1.0), 'm-',label="d = 20")

pl.plot(v_vals, finalV(v_vals, d=0.5, lam=1.0, kind="approx"), 'b--')
pl.plot(v_vals, finalV(v_vals, d=2.0, lam=1.0, kind="approx"), 'g--')
pl.plot(v_vals, finalV(v_vals, d=5.0, lam=1.0, kind="approx"), 'r--')
pl.plot(v_vals, finalV(v_vals, d=20.0, lam=1.0, kind="approx"), 'm--')

pl.plot([0,800], [0, 800], 'k:')

pl.xlabel(r"$v_i$ [km/s]")
pl.ylabel(r"$v_f$ [km/s]")

pl.legend(loc="best")

pl.show()




print " Calculating for Germanium..."
sig_corr_Fe = v_vals*0.0
ER_corr_Fe = v_vals*0.0
for i, v in enumerate(v_vals):
    sig_corr_Fe[i] = quad(lambda x: DMU.calcSIFormFactor(x*(1e6/(3e5*3e5))*ERmax(m_x, 0.9315*73, v), 73), 0, 1)[0]
    ER_corr_Fe[i] = quad(lambda x: 2.0*x*DMU.calcSIFormFactor(x*(1e6/(3e5*3e5))*ERmax(m_x, 0.9315*73, v), 73), 0, 1)[0]

print " Calculating for Oxygen..."
sig_corr_O = v_vals*0.0
ER_corr_O = v_vals*0.0
for i, v in enumerate(v_vals):
    sig_corr_O[i] = quad(lambda x: DMU.calcSIFormFactor(x*(1e6/(3e5*3e5))*ERmax(m_x, 0.9315*16, v), 16), 0, 1)[0]
    ER_corr_O[i] = quad(lambda x: 2.0*x*DMU.calcSIFormFactor(x*(1e6/(3e5*3e5))*ERmax(m_x, 0.9315*16, v), 16), 0, 1)[0]

f, axarr = pl.subplots(2)


ax1 = axarr[0]
ax2 = axarr[1]
#Subplot 1 - cross section correction
ax1.plot(v_vals,sig_corr_O, color='blue',label='O', linewidth=1.5)
ax1.plot(v_vals,sig_corr_Fe, color='green', label='Fe', linewidth=1.5)

#ax1.plot(v_vals, 1-(v_vals/(800*2.3))**2, 'b--')


#ax1.set_xticks([])
ax1.set_ylabel('Cross section correction')
#ax1.set_ylim(0, 1.1)
#ax1.axhline(1.0, color='k', linestyle=':')

#Subplot 2 - recoil energy correction
ax2.plot(v_vals,ER_corr_O, color='blue',label='O', linewidth=1.5)
ax2.plot(v_vals,ER_corr_Fe, color='green', label='Pb', linewidth=1.5)

#ax2.plot(v_vals, 1-(v_vals/(800*2.1))**2, 'b--')
#ax2.plot(v_vals, 1-(v_vals/(800*0.5))**2 , 'g--')
#ax2.plot(v_vals**2, np.exp(-(v_vals/(800*0.40))**2) + (v_vals/(8000*0.40))**2*np.exp(-(v_vals/(800))**2) , 'g--')
#ax2.plot(v_vals, np.exp(-(v_vals/(800*0.5))**2) , 'g--')

ax2.set_xlabel('v [km/s]', fontsize=20.0)
ax2.set_ylabel(r'$C_i(v) = \sigma_i(v)/\sigma_i(0)$', fontsize=22.0)
ax2.set_ylim(0, 1.1)
ax2.axhline(1.0, color='k', linestyle=':')
ax2.legend(loc='best', frameon=False)

#f.subplots_adjust(hspace=0.1)
pl.savefig("test.pdf", bbox_inches="tight")

n_Ge = (5.3/73.0)*6.02e23

pl.figure()
pl.plot(v_vals,1*n_Ge*1e6*1e-28*(73**4)*(73*0.9315)*(v_vals/3e5)**2*ER_corr_Fe[i])

pl.show()
