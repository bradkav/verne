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

import matplotlib.pyplot as plt

import verne


# Some random functions...
def q(A, E_R):
    mN = 1e3*A*0.9315 #in MeV
    
    return np.sqrt(2*mN*E_R*1e-3) #in MeV

def mu(mX, mA):
    return mX*mA*1.0/(mX+mA)

def ERmax(mX, mA, v):
    return 2*(mu(mX, mA)*v)**2/mA


m_x = 1e5

v_vals = np.linspace(0.1, 800, 100)

#Calculating correction factors...

verne.loadFFcorrections(m_x, interaction = "millicharge")
#verne.loadFFcorrections(m_x, interaction = "SI")

plt.figure()

for key in verne.isoID:
    if ("_A" not in key):
        iso = verne.isoID[key]
    
        #print(key, iso)
        C_vals = verne.corr_interp[iso](v_vals)
        #print(C_vals)
        plt.plot(v_vals, C_vals, label=key, color = "C" + str(iso))
        
        m_A = 0.93*verne.Avals[iso]
        mu_A = m_A*m_x/(m_A + m_x)
        q_max = 2*mu_A*v_vals/3e5
        x_s = (verne.q_screen/q_max)**2
        
        v_min = 3e5*verne.q_screen/(2*mu_A)
        #print(v_min)
        
        #plt.plot(v_vals, np.log(1/x_s), color = "C" + str(iso), linestyle='--')
        
        plt.axvline(v_min, linestyle=':', color='k')
        
    

plt.xlabel(r'$v$ [km/s]', fontsize=20.0)
plt.ylabel(r'$C_i(v)$', fontsize=12.0)
plt.ylim(0, 1.1)
plt.axhline(1.0, color='k', linestyle=':')
plt.legend(loc='best')

plt.show()