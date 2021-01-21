import numpy as np
from numpy import pi
from scipy.integrate import simps,quad, cumtrapz
import verne
from LabFuncs import *
import utils
from scipy.special import erf
from scipy.interpolate import interp1d
import MaxwellBoltzmann as MB
import argparse
import os.path
from tqdm import tqdm

from WIMpy import DMUtils as DMU_standard
import DMUtils_migdal as DMU

from matplotlib import pyplot as plt

try:
    from tqdm import tqdm
except:
    def tqdm(x):
        return x


#Parse the arguments!
parser = argparse.ArgumentParser(description='...')
parser.add_argument('-m_x','--m_x', help='DM mass in GeV', type=float,default = 1e5)
parser.add_argument('-sigma_p','--sigma_p', help='DM-nucleon cross section, sigma_p in cm^2', type=float, required=True)
#parser.add_argument('-latitude','--latitude', help='Latitude (in degerees) of detector', type=float, required=True)
parser.add_argument('-interaction', '--interaction', help='Interaction type: `SI` or `SD`', type=str, default="SI")

args = parser.parse_args()
m_x = args.m_x
sigma_p = args.sigma_p
interaction = args.interaction

lat = 45.2 #degrees

A_Ge = 72.64

t_exp = 30.0
#Get Julian date of exposure
#JulianDay(month, day, year, hour)
t0 = JulianDay(3, 15, 2020, 0)
t1 = 0.0 #Start time is 14 minutes past 22hrs
t2 = t1 + t_exp #Total run time is 5.31 hours

gamma_grid = np.pi*np.linspace(0, 1, 11)


#Average over 30 days

Nvals = 10001
tvals = t0 + np.linspace(t1, t2, Nvals)

gammavals = np.zeros(Nvals)

#Calculate gamma from the LabVelocity
for i in range(Nvals):
    vs = -LabVelocity(tvals[i], lat, 0)
    vs_hat = vs/np.sqrt(np.sum(vs**2))
    rdet_hat = np.asarray([0,0,1.0])
    
    gammavals[i] = np.arccos(np.dot(vs_hat, rdet_hat))


#Load velocity distribution from file
#def getVelDist(m_x, sigma_p, gamma_ind):
#    Ngamvals = 11
#    Nvvals = 62
#    
#    fname = "../results/veldists/f_" + interaction + "_MOD_mx" + '{0:.3f}'.format(m_x) + "_lsig" + '{0:.2f}'.format(np.log10(sigma_p)) + ".txt"
#    
#    gamma_vals1, vvals1, fvals1 = np.loadtxt(fname, unpack=True)
#    vvals = vvals1[gamma_ind*Nvvals:(gamma_ind+1)*Nvvals]
#    fvals = fvals1[gamma_ind*Nvvals:(gamma_ind+1)*Nvvals]
#    return vvals, fvals
    

fname = "results/veldists/f_" + interaction + "_MOD_mx" + '{0:.3f}'.format(m_x) + "_lsig" + '{0:.2f}'.format(np.log10(sigma_p)) + ".txt"
gamma_vals1, vvals1, fvals1 = np.loadtxt(fname, unpack=True)
v_grid = np.reshape(vvals1, (11, 62))
f_grid = np.reshape(fvals1, (11, 62))
   
v_cut = 20.0
   
v_list = np.linspace(0.0, 775.0, 1000)
f_grid_new = np.zeros((11, len(v_list)))
   
for i in range(11):
    inds = np.argsort(v_grid[i,:])
    v1 = v_grid[i, inds]
    f1 = f_grid[i, inds]
    #print(i, np.trapz(f1,v1))
    f_grid_new[i,:] = np.interp(v_list, v1, f1, left=0, right=0)

f_avg = 0.0*v_list

for j, v in enumerate(v_list):
    f_vals = np.interp(gammavals, gamma_grid, f_grid_new[:, j], left=0, right=0)
    f_avg[j] = np.trapz(f_vals, tvals)/t_exp
#DON'T CORRECT THE NORMALISATION HERE!!!

f_clip = 0.0*f_avg
f_clip[v_list > v_cut] = f_avg[v_list > v_cut]
#v_list, f_avg
#print(np.trapz(f_avg, v_list))
    
#plt.figure()
#
#for i in range(11):
#    plt.plot(v_grid[i,:], f_grid[i,:])
#
#plt.plot(v_list, f_avg, color='k', lw=2)
#plt.plot(v_list, f_clip, color='k', linestyle=':', lw=2)
#

#print(f_avg)

N_ER = 200
N_EEM = 200

ER_list = np.geomspace(1e-5, 1e0, N_ER)
EEM_list = np.geomspace(1e-5, 1e0, N_EEM)

ER_grid, EEM_grid = np.meshgrid(ER_list, EEM_list, indexing='ij')

#--------------------------------------------------
#--------- Standard Elastic Nuclear Recoils -------
#--------------------------------------------------

eta_list = np.trapz(f_clip/(v_list+1e-30), v_list) - cumtrapz(f_clip/(v_list+1e-30), v_list, initial = 0)
eta_func = interp1d(v_list, eta_list, bounds_error=False, fill_value=0.0)

#plt.figure()
#
#plt.plot(v_list, eta_list)


R_NR_list = DMU_standard.dRdE_generic(ER_list, A_Ge, 0, m_x, sigma_p, eta_func)
R_NR_list_noscatter = DMU_standard.dRdE_standard(ER_list, A_Ge, 0, m_x, sigma_p)

outfile = "results/spectra/Elastic/Elastic_" + interaction + "_MOD_mx" + '{0:.3f}'.format(m_x) + "_lsig" + '{0:.2f}'.format(np.log10(sigma_p)) + ".txt"
headertxt = "Elastic nuclear recoil spectrum for " + interaction + " interactions. Ge detector at Modane (d=1700m, latitude=45.2deg)"
headertxt += "\nColumns: Nuclear recoil energy [keV], dR/dE (w/ Earth-shielding) [events/keV/kg/day]"   
#Output format - short
#np.savetxt(outfile, list(zip(ER_list, R_NR_list, R_NR_list_noscatter)), header=headertxt, fmt='%.4e')
np.savetxt(outfile, list(zip(ER_list, R_NR_list)), header=headertxt, fmt='%.3e')
   
#plt.figure()
#
#plt.loglog(ER_list, R_NR_list_noscatter, '--')
#plt.loglog(ER_list, R_NR_list)
#
#plt.show()

#--------------------------------
#--------- Migdal Recoils -------
#--------------------------------
   
f_func = interp1d(v_list, f_clip, bounds_error=False, fill_value=0.0)
   
#integ = lambda x: DMU.d2RdERdEe(x, Etot-x, m_x, sigma_p, A_Ge, n, speeddist=f_func, SD=nucleon, force_quad=False)
   
d2R = ER_grid*0.0

#Only run if the velocity distribution is non-zero
if (np.sum(f_clip) > 1e-40):
    
    for i in tqdm(range(N_ER)):
        for j in range(N_EEM):
            d2R[i,j] = DMU.d2RdERdEe(ER_list[i], EEM_list[j], m_x, sigma_p, A_Ge, n_sep=3, speeddist=f_func, force_quad=False)
else:
    print(">Rate is zero...")
    
out_data = list(zip(ER_grid.flatten(), EEM_grid.flatten(), d2R.flatten()))

htext = "Migdal Spectrum (with Earth-shielding), including only n = 3 shell electrons. Ge detector at Modane (d=1700m, latitude=45.2deg)."
htext += "\nColumns: E_NR [keV], E_e [keV], d2R/dE_NR/dE_e [evts/kg/day/keV^2]"

np.savetxt("results/spectra/Migdal2D/Migdal2D_" +  interaction + "_MOD_mx" + '{0:.3f}'.format(m_x) + "_lsig" + '{0:.2f}'.format(np.log10(sigma_p)) + ".txt",out_data, header=htext,  fmt='%.3e')

#1D spectrum

dRdEe = np.trapz(d2R, ER_list, axis = 0)

out_data = list(zip(EEM_list, dRdEe))

htext = "Migdal Spectrum, including only n = 3 shell electrons. Integrated over nuclear recoil energies. Ge detector at Modane (d=1700m, latitude=45.2deg)."
htext += "\nColumns: E_e [keV], dR/dE_e [evts/kg/day/keV]"

np.savetxt("results/spectra/Migdal1D/Migdal1D_" +  interaction + "_MOD_mx" + '{0:.3f}'.format(m_x) + "_lsig" + '{0:.2f}'.format(np.log10(sigma_p)) + ".txt",out_data, header=htext,  fmt='%.3e')

#plt.figure()
#
#plt.loglog(EEM_list, dRdEe)
#
#plt.xlim(2e-3, 1.1)
#
#plt.savefig("/Users/bradkav/Desktop/Migdal.pdf")
#
#plt.show()
#
#----------------------------


#Calculate velocity integral from an interpolation function
#defining f(v) and a maximum speed vmax
#def calcEta_final(v, interpfun, vmax):
#    return quad(lambda x: interpfun(x)/x, v, vmax*1.1)[0]

