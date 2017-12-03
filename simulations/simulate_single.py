from verne import *
from tqdm import tqdm
import WIMpy.DMUtils as DMU

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

R_E = 6371.0 #km

v0 = np.asarray([0, 0, 750.0])
x0 = np.asarray([0,0,0])
Nsteps = 1000000

v_list = np.zeros((3, Nsteps))
vmag_list = np.zeros(Nsteps)

x_list = np.zeros((3, Nsteps))

dt = 1e-11

tvals = np.arange(Nsteps)*dt

rvals = np.zeros(1000)

l = (m_N/m_x)/Lambda

print "lambda = ",l

svals = (2.0/3)*l*(m_N/m_x)*0.5*tvals*v0[2]**3*(2+tvals*v0[2]*l)/(1 + tvals*v0[2]*l)**2
svals2 = (2.0/3)*l*(m_N/m_x)*0.5*tvals*v0[2]**3*(2+tvals*v0[2]*l)/(1 + tvals*v0[2]*l)**3

#svals2 = (2.0/3)*l*(m_N/m_x)*tvals*v0[2]**3*(3 + tvals*v0[2]*l*(3+tvals*l*v0[2]))/(3.0*(1+tvals*v0[2]*l)**3)

def gaussian(x, mu, sig):
    return (np.sqrt(2*np.pi*sig**2)**-1)*np.exp(-np.power(x - mu, 2.) / (2 * np.power(sig, 2.)))

v_current = v0
x_current = x0
total_distance = 0
for i in tqdm(range(Nsteps)):
    v_list[:,i] = v_current
    x_list[:,i] = x_current
    vmag_list[i] = np.sqrt(np.sum(v_current**2))
    x_current, v_current, deltax = step(x_current, v_current, dt)
    total_distance += deltax

print " Total z-distance travelled: ", x_list[2, -1]
print " Total distance travelled:   ", total_distance
pl.figure()
pl.plot(x_list[2,:], vmag_list**2/v0[2]**2, label="MC")
pl.plot(x_list[2,:], np.exp(-2*l*x_list[2,:]))
#pl.plot(tvals, v0[2]/(1+v0[2]*l*tvals), label="Analytic (approx)")
pl.legend(loc="best")
pl.ylabel(r'$E/E_0$')
pl.xlabel('d [km]')

pl.figure()
pl.plot(x_list[2,:], vmag_list, label="MC")
pl.plot(x_list[2,:], v0[2]/(1+v0[2]*l*tvals), label="Analytic (approx)")
pl.legend(loc="best")
pl.ylabel(r'$|v| [km/s]$')
pl.xlabel('d [km]')

pl.figure()
pl.plot(tvals, np.sqrt(v_list[0,:]**2 + v_list[1,:]**2),label=r"$v_T$")
pl.plot(tvals, np.sqrt((1.0/3)*l*(m_N/m_x)*tvals*v0[2]**3*(2+tvals*v0[2]*l)/(1 + tvals*v0[2]*l)**2), label="Analytic")
pl.plot(tvals, np.sqrt((1.0/3)*l*(m_N/m_x)*tvals*v0[2]**3*(2+tvals*v0[2]*l)/(1 + tvals*v0[2]*l)**3), label="Analytic (fudged)")
#pl.semilogy(tvals, np.abs(v_list[1,:]),label=r"$v_y$")
#pl.semilogy(tvals, np.abs(v_list[2,:]),label=r"$v_z$")

pl.ylabel(r'$v_i [km/s]$')
pl.xlabel('t [s]')
pl.legend(loc="best")

pl.figure()
pl.semilogy(tvals, np.sqrt(x_list[0,:]**2 + x_list[1,:]**2),label=r"$r_T$")
#pl.semilogy(tvals, np.abs(x_list[1,:]),label=r"$y$")
pl.semilogy(tvals, np.abs(x_list[2,:]),label=r"$z$")
pl.ylabel(r'$r_i [km]$')
pl.xlabel('t [s]')
pl.legend(loc="best")

pl.axhline(2*R_E, linestyle='--')

pl.ylim(1e-5,1e5)

pl.figure()
pl.plot(np.sqrt(x_list[0,:]**2 + x_list[1,:]**2),x_list[2,:])
pl.xlabel('r_T [km]')
pl.ylabel('z [km]')
#pl.xlim(0, 2.5)
#pl.ylim(0, 2.5)

print " ", 1e3*np.sqrt(x_list[0,-1]**2 + x_list[1,-1]**2), " metres of deflection"

print " ", np.sqrt(x_list[0,-1]**2 + x_list[1,-1]**2)/x_list[2,-1], " radians deviation"

pl.show()
    