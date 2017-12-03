import numpy as np
import WIMpy.DMUtils as DMU
from numpy import pi
from scipy.integrate import quad, dblquad
import verne_analytic as VA
from verne_analytic import isoID
from skmonaco import mcquad
from timeit import default_timer as timer
import sys
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


if (len(sys.argv) != 2):
    print " CalculateVT.py requires 1 argument - e.g. CalculateVT.py LOG10SIG"
    print " Exiting..."
    sys.exit()
    
#print " I've changed the grid system!"
#sys.exit()
    
lsigstr = sys.argv[1]
sigma_p = 10**float(lsigstr)

m_x = 1e5


VA.loadIsotopes()
VA.loadPhiInterp()

Nvals = 50

Nthetavals = 10

thetalist = np.linspace(np.pi/2.0, np.pi, Nthetavals)

#theta = np.pi*(0.5+0.5)
depth = 10.6
#sigma_p = 1e-28
#m_x = 1e5

v_initial = np.linspace(0, 800, Nvals)
v_final_lab = np.zeros((Nvals, Nthetavals))
v_final_full = np.zeros((Nvals, Nthetavals))
calcVf = np.vectorize(VA.calcVfinal_full)


for j, theta in enumerate(tqdm(thetalist)):
    for i, v in enumerate(v_initial):
        #v_final[i] = VA.calcVfinal_full(v, theta,  depth, sigma_p, m_x, target="atmos")
        v_final_lab[i,j] = VA.calcVfinal_full(v, theta,  depth, sigma_p, m_x, target="full")
        v_final_full[i,j] = VA.calcVfinal_shield(v_final_lab[i, j], theta,  depth, sigma_p, m_x)

vgrid, tgrid = np.meshgrid(v_initial, thetalist)

output = zip(vgrid.flatten(), tgrid.flatten(), np.clip((v_final_full.T).flatten(),0,800.0))

np.savetxt("results/VT/VT_lsig=" + lsigstr + ".txt", output)
