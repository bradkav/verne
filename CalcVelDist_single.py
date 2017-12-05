import numpy as np
from scipy.interpolate import interp1d
import verne
import MaxwellBoltzmann as MB
import argparse

from matplotlib import pylab as pl

#Parse the arguments!
parser = argparse.ArgumentParser(description='...')
parser.add_argument('-m_x','--m_x', help='DM mass in GeV', type=float,default = 1e5)
parser.add_argument('-sigma_p','--sigma_p', help='DM-nucleon cross section, sigma_p in cm^2', type=float, required=True)
parser.add_argument('-gamma','--gamma_by_pi', help='Average incoming DM angle, gamma, ***IN UNITS OF PI***. Gamma = pi corresponds to DM particles coming mostly from overhead.', type=float, required=True)
parser.add_argument('-target','--target', help='Target to propagate through. `earth`, `atmos`, `shield` or `full`.', type=str, default="full")
args = parser.parse_args()
m_x = args.m_x
sigma_p = args.sigma_p
gamma_by_pi = args.gamma_by_pi
gamma = np.pi*gamma_by_pi
target = args.target

#Stanford lab
depth = 10.6 #metres

print "   "
print "    Calculating for..."
print "        m_x/GeV:", m_x
print "        sigma_p/cm^2:", sigma_p
print "        gamma/pi:", gamma_by_pi
print "        target:", target
print " "

#Initialise verne
verne.loadIsotopes()
verne.loadFFcorrections(m_x)

#Calculate the maximum initial speed as a function of incoming angle theta
Nvals = 1001
thetavals = np.linspace(0, np.pi, Nvals)

v_e = np.sqrt(2.0)*MB.sigmav
vesc = MB.vesc

a = 1.0
b = 2*v_e*(-np.sin(gamma)*np.sin(np.pi-thetavals) + np.cos(gamma)*np.cos(np.pi-thetavals))
c = v_e**2 - vesc**2

v_initial_max = (-b + np.sqrt(b**2 - 4*a*c))/(2.0*a)

#Calculate the maximum final speed as a function of incoming angle theta
v_final_max = 0.0*v_initial_max
for i in range(Nvals):
    v_final_max[i] = verne.calcVfinal_full(v_initial_max[i], thetavals[i],  depth, sigma_p, m_x, target)

pl.figure()
pl.plot(thetavals/np.pi,v_final_max)
pl.show()

#Calculate interpolation function for max final speed
vmax = np.max(v_final_max)
vfinal_interp = interp1d(thetavals, v_final_max, kind='linear', bounds_error=False, fill_value=0)


#Tabulate values of speed distribution
v_th = 1.0 #Lowest speed to consider (don't go lower than 1 km/s, other the calculation of derivatives is messed up...)

#Generate a list of sampling values for v (with some very close to v_th)
vlist = np.logspace(np.log10(v_th), np.log10(0.25*vmax), 20)
vlist = np.append(vlist, np.linspace(0.15*vmax, 0.99*vmax, 40))
#vlist = np.linspace(v_th, 0.999*vmax, 50)
vlist = np.sort(vlist)
f_final = 0.0*vlist
for i in range(len(vlist)):
    f_final[i] = verne.CalcF(vlist[i], gamma, depth, sigma_p, m_x, target, vfinal_interp)

#Add on the final point
vlist = np.append(vlist, vmax)
f_final = np.append(f_final, 0.0)

pl.figure()
pl.plot(vlist, f_final, "+-")

pl.figure()
pl.semilogy(vlist, f_final, "+-")

pl.show()

#Output to file
fname = "results/veldists/f_lmx" + '{0:.1f}'.format(np.log10(m_x)) + "_lsig" + '{0:.1f}'.format(np.log10(sigma_p)) + "_gam" + '{0:.1f}'.format(gamma_by_pi) + ".txt"
headertxt = "mx [GeV]: " + str(m_x) + "\nsigma_p [cm^2]: " + str(sigma_p) + "\ngamma/pi: " + str(gamma_by_pi) + "\ndepth [m]: " + str(depth) + "\ntarget: " + target
headertxt += "\nColumns: v [km/s], f(v) [s/km]"

np.savetxt(fname, np.transpose([vlist, f_final]), header=headertxt)
    