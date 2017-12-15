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

ax1.fill(alt_data[:,0], alt_data[:,1], color='DarkOrange', alpha=0.25)
ax1.fill(nucleus_data[:,0], nucleus_data[:,1], color='grey', alpha=0.25)

ax1.plot(upper_old[:,0], upper_old[:,1], 'k-', linewidth=2.0)

ax1.plot(upper_new[0:-1,0], upper_new[0:-1,1], 'b--',linewidth=2.0)
ax1.plot(upper_new[-2:,0], upper_new[-2:,1], 'b-',linewidth=2.0)
ax1.plot(upper_new_nucleus[0:-1,0], upper_new_nucleus[0:-1,1], 'r--',linewidth=2.0)
ax1.plot(upper_new_nucleus[-2:,0], upper_new_nucleus[-2:,1], 'r-',linewidth=2.0)

ax1.set_yscale("log")
ax1.set_xscale("log")


ax1.fill_between(mvals, np.vectorize(calc_lower)(mvals), upper_old_interp(mvals), alpha=0.25, color='grey',edgecolor='black')
ax1.fill_between(mvals2, upper_old_interp(mvals2), upper_new_interp(mvals2), alpha=0.4, color='DarkBlue',edgecolor='black')
ax1.fill_between(mvals5, upper_old_interp(mvals5), 10**upper_new_nucleus_interp(np.log10(mvals5)), alpha=0.4, color='DarkRed',edgecolor='black')
ax1.set_ylim(1e-47, 1e-17)

ax1.set_xlabel(r"DM mass $m_\chi \,\,[\mathrm{GeV}]$")
ax1.set_ylabel(r"DM-nucleon cross section $\sigma_p^\mathrm{SI}\,\,[\mathrm{cm}^2]$")
ax1.xaxis.set_minor_locator(LogLocator(base=10.0, subs=(10,)))
ax1.yaxis.set_minor_locator(LogLocator(base=10.0, subs=(10,100,1000)))

txtfont = 12.0

ax1.text(10, 1e-19, "High-altitude", color='DarkOrange',fontsize=txtfont)
ax1.text(1e6, 1e-35, "Direct detection",fontsize=txtfont)
ax1.text(1e4, 1e-45, "Xenon1T",fontsize=txtfont)
ax1.text(1e0, 1e-38, "CRESST-III",fontsize=txtfont)
#ax1.text(5e-2, 1e-26, "CRESST 2017\nsurface",fontsize=txtfont)
ax1.text(1e6, 1e-21, "CDMS-I (this work)", color='DarkBlue',fontsize=txtfont)
ax1.text(0.5e1, 0.3e-25, "CRESST 2017 surface\n(this work)", color='Red',fontsize=txtfont)


pl.savefig('../plots/Constraints1.pdf', bbox_inches='tight',fontsize=txtfont)
pl.show()





