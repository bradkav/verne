import numpy as np
from scipy.interpolate import interp1d
from matplotlib.ticker import MultipleLocator, FormatStrFormatter,LogLocator
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

alt_data = np.loadtxt("../results/constraints/high-altitude.txt")

IC_data = np.loadtxt("../results/constraints/IceCube.txt")

CR_data = np.loadtxt("../results/constraints/CosmicRays.txt")

DD_data = np.loadtxt("../results/constraints/directdetection_all.txt")

EarthHeat_data = np.loadtxt("../results/constraints/EarthHeat.txt")

CRESST_data = np.loadtxt("../results/constraints/CRESST-III.txt")
CRESST_data[:,1] *= 1e-38

XENON_data = np.loadtxt("../results/constraints/x1t_firstdata_result_edit.txt", skiprows=1, usecols=(0,1), delimiter=',')
XENON_data[:,1] *= 1e-45

nucleus_data = np.loadtxt("../results/constraints/nucleus.txt")

upper_old = np.loadtxt("../results/constraints/upper_edge.txt")
upper_old_interp = interp1d(upper_old[:,0], upper_old[:,1], fill_value="extrapolate", kind="linear")


upper_new = np.loadtxt("../results/constraints/CDMS_thiswork.txt")
upper_new_interp = interp1d(upper_new[:,0], upper_new[:,1], fill_value="extrapolate")

upper_new_nucleus = np.loadtxt("../results/constraints/nucleus_thiswork.txt")
upper_new_nucleus_interp = interp1d(np.log10(upper_new_nucleus[:,0]), np.log10(upper_new_nucleus[:,1]), kind="linear")

CRESST_interp = interp1d(CRESST_data[:,0], CRESST_data[:,1], bounds_error=False,fill_value=1e10)
XENON_interp = interp1d(XENON_data[:,0], XENON_data[:,1]*1e45, bounds_error=False, fill_value=1e10)


mvals = np.logspace(np.log10(0.36), np.log10(3e16), 100)
mvals2 = np.logspace(np.log10(1e1), np.log10(1e15), 100)
mvals4 = np.logspace(np.log10(1e0), np.log10(1e5), 20)
mvals5 = np.logspace(np.log10(1e0), np.log10(1e8), 20)

def calc_lower(m):
    return np.minimum(CRESST_interp(m), XENON_interp(m)*1e-45)

pl.figure()
ax1 = pl.gca()

ax1.fill(DD_data[:,0], DD_data[:,1], facecolor='Grey', alpha=0.5, linestyle='-', edgecolor='k',linewidth=2.0, zorder=0)
ax1.fill(alt_data[:,0], alt_data[:,1], facecolor='Orange', alpha=0.25)
ax1.add_patch(mpl.patches.Polygon(alt_data, facecolor='None', edgecolor='k', linestyle='--'))
ax1.fill(IC_data[:,0], IC_data[:,1], facecolor='Red', alpha=0.25)
ax1.add_patch(mpl.patches.Polygon(IC_data, facecolor='None', edgecolor='k', linestyle=':'))
ax1.fill(CR_data[:,0], CR_data[:,1], facecolor='Green', alpha=0.25)
ax1.add_patch(mpl.patches.Polygon(CR_data, facecolor='None', edgecolor='k', linestyle='-'))
ax1.fill(EarthHeat_data[:,0], EarthHeat_data[:,1], facecolor='Cyan', alpha=0.25)
ax1.add_patch(mpl.patches.Polygon(EarthHeat_data, facecolor='None', edgecolor='k', linestyle='-.'))


ax1.set_yscale("log")
ax1.set_xscale("log")


ax1.set_ylim(1e-47, 1e-17)
ax1.set_xlim(1e-2, 1e18)

ax1.set_xlabel(r"DM mass $m_\chi \,\,[\mathrm{GeV}]$")
ax1.set_ylabel(r"DM-nucleon cross section $\sigma_p^\mathrm{SI}\,\,[\mathrm{cm}^2]$")
ax1.xaxis.set_minor_locator(LogLocator(base=10.0, subs=(10,)))
ax1.yaxis.set_minor_locator(LogLocator(base=10.0, subs=(10,100,1000)))

txtfont = 12.0

ax1.text(1.5e10, 1e-19, "High-altitude", color='DarkOrange',fontsize=txtfont)
ax1.text(1e6, 1e-35, "Direct detection",fontsize=txtfont)
ax1.text(1, 1e-19, "Cosmic Rays", color='DarkGreen',fontsize=txtfont)
ax1.text(10, 5e-31, "Earth Heat", color='DarkCyan',fontsize=txtfont)
ax1.text(1e7, 5e-25, "IceCube", color='DarkRed',fontsize=txtfont)

pl.savefig('../plots/Constraints2.pdf', bbox_inches='tight',fontsize=txtfont)
pl.show()





