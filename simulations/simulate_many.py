from verne import *
from tqdm import tqdm

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
Nsteps = 1000

v_list = np.zeros((3, Nsteps))
vmag_list = np.zeros(Nsteps)

x_list = np.zeros((3, Nsteps))

dt = 0.001

tvals = np.arange(Nsteps)*dt

rvals = np.zeros(1000)

l = (m_N/m_x)/Lambda

svals = (2.0/3)*l*(m_N/m_x)*0.5*tvals*v0[2]**3*(2+tvals*v0[2]*l)/(1 + tvals*v0[2]*l)**2
svals2 = (2.0/3)*l*(m_N/m_x)*tvals*v0[2]**3/(1 + tvals*v0[2]*l)

def gaussian(x, mu, sig):
    return (np.sqrt(2*np.pi*sig**2)**-1)*np.exp(-np.power(x - mu, 2.) / (2 * np.power(sig, 2.)))

for k in range(1000):
    v_current = v0
    x_current = x0
    for i in tqdm(range(Nsteps)):
        v_list[:,i] = v_current
        x_list[:,i] = x_current
        vmag_list[i] = np.sqrt(np.sum(v_current**2))
        x_current, v_current = step(x_current, v_current, dt)

    rvals[k] = v_list[0,-1]**2 + v_list[1,-1]**2
    rvals[k] = v_list[0,-1]

    pl.show()
    
    
print np.mean(rvals**2)
print 0.5*svals[-1]
print 0.5*svals2[-1]

r_ext = np.maximum(np.abs(np.min(rvals)), np.abs(np.max(rvals)))
binvals = np.linspace(-r_ext, r_ext, 101)

pl.figure()
pl.hist(rvals, binvals)
#pl.hist(rvals, 50)
pl.plot(binvals, 1000*(r_ext/50.0)*(gaussian(binvals, 0, np.sqrt(0.5*svals[-1]))), linewidth=2.0)
pl.show()


vmag_list = vmag_list[::10]
v_list = v_list[:,::10]
x_list = x_list[:,::10]
tvals = tvals[::10]
    
l = (m_N/m_x)/Lambda
    
pl.figure()
pl.plot(tvals, vmag_list, label="MC")
pl.plot(tvals, v0[2]/(1+v0[2]*l*tvals), label="Analytic (approx)")
pl.legend(loc="best")
pl.ylabel(r'$|v| [km/s]$')
pl.xlabel('t [s]')

pl.figure()
pl.plot(tvals, np.sqrt(v_list[0,:]**2 + v_list[1,:]**2),label=r"$v_T$")
pl.plot(tvals, np.sqrt((1.0/3)*l*(m_N/m_x)*tvals*v0[2]**3*(2+tvals*v0[2]*l)/(1 + tvals*v0[2]*l)**2), label="Analytic (approx)")
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

print np.sqrt(x_list[0,Nsteps/20]**2 + x_list[1,Nsteps/20]**2)

pl.show()
    
    