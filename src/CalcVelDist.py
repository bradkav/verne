import numpy as np
from scipy.interpolate import interp1d

import MaxwellBoltzmann as MB
import argparse

try:
    from tqdm import tqdm
except:
    def tqdm(x):
        return x
    
results_dir = "results/"

#Parse the arguments!
parser = argparse.ArgumentParser(description='...')
parser.add_argument('-m_x','--m_x', help='DM mass in GeV', type=float,default = 1e5)
parser.add_argument('-sigma_p','--sigma_p', help='DM-nucleon cross section, sigma_p in cm^2', type=float, required=True)
parser.add_argument('-loc','--location', help='Detector location to consider.', type=str, required=True)
parser.add_argument('-interaction', '--interaction', help='Interaction type: `SI` or `SD`', type=str, default="SI")
args = parser.parse_args()
m_x = args.m_x
sigma_p = args.sigma_p
loc = args.location
interaction = args.interaction

if (interaction == "SI"):
    import verne as verne
elif (interaction == "SD"):
    import verne_SD as verne
else:
    print(">Unknown interaction type <", interaction, ">...")
    exit()

#exit()


if (loc == "SUF"):
    depth = 10.6 #metres
elif (loc == "MPI"):
    depth = 0.3 #metres
elif (loc == "EDE"):
    depth = 0.0 #metres
elif (loc == "surface"):
    depth = 0.0 #metres
elif (loc == "MOD"):
    depth = 1700.0

target = loc


print(" ")
print(">Calculating for...")
print(">    m_x/GeV:", m_x)
print(">    sigma_p/cm^2:", sigma_p)
#print "        gamma/pi:", gamma_by_pi
print(">    detector at :", loc)
print(" ")

#Initialise verne
verne.loadIsotopes()
verne.loadFFcorrections(m_x)

#Calculate the maximum initial speed as a function of incoming angle theta
Nvals = 501
thetavals = np.linspace(0, np.pi, Nvals)

v_e = np.sqrt(2.0)*MB.sigmav
vesc = MB.vesc

v_th = 1e0 #Lowest speed to consider (don't go lower than v_th km/s, other the calculation of derivatives is messed up...)  

def getVelDist(gamma):
    
    print(">    Calculating maximum final speed...")
    a = 1.0
    b = 2*v_e*(-np.sin(gamma)*np.sin(np.pi-thetavals) + np.cos(gamma)*np.cos(np.pi-thetavals))
    c = v_e**2 - vesc**2
    v_initial_max = (-b + np.sqrt(b**2 - 4*a*c))/(2.0*a)

    #Calculate the maximum final speed as a function of incoming angle theta
    v_final_max = 0.0*v_initial_max
    for i in tqdm(range(Nvals)):
        v_final_max[i] = verne.calcVfinal_full(v_initial_max[i], thetavals[i],  depth, sigma_p, m_x, target)

    #Calculate interpolation function for max final speed
    vmax = np.max(v_final_max)
    if (vmax < v_th):
        return np.linspace(0, 1, Nv), np.zeros(Nv)
    
    vfinal_interp = interp1d(thetavals, v_final_max, kind='linear', bounds_error=False, fill_value=0)
    
    print(">    Calculating final speed distribution...")

    #Tabulate values of speed distribution

    #Generate a list of sampling values for v (with some very close to v_th)
    vlist = np.logspace(np.log10(v_th), np.log10(0.25*vmax), Nv1)    #20
    vlist = np.append(vlist, np.linspace(0.15*vmax, 0.999*vmax, Nv2)) #40
    vlist = np.append(vlist, 19.9)
    #vlist = np.linspace(v_th, 0.999*vmax, 50)
    vlist = np.sort(vlist)
    f_final = 0.0*vlist
    for i in tqdm(range(len(vlist))):
        f_final[i] = verne.CalcF(vlist[i], gamma, depth, sigma_p, m_x, target, vfinal_interp)

    #Add on the final point
    vlist = np.append(vlist, vmax)
    f_final = np.append(f_final, 0.0)

    return vlist, f_final
    
    
#Loop over gamma values
N_gamma = 11 #Do 21 in future
Nv1 = 20 #Near the velocity threshold
Nv2 = 40 #Everywhere else
Nv = Nv1 + Nv2 + 1 + 1 #Add an extra one for 20 km/s
gamma_list = np.linspace(0, 1.0, N_gamma)
vgrid = np.zeros((N_gamma, Nv))
fgrid = np.zeros((N_gamma, Nv))
for j in range(N_gamma):
    print(">Calculating for gamma/pi = ", gamma_list[j],"...")
    vgrid[j,:], fgrid[j,:] = getVelDist(gamma_list[j]*np.pi)

gamma_rep = np.repeat(gamma_list, Nv)

#Output to file
fname = results_dir + "veldists/f_" + interaction + "_" + loc + "_mx" + '{0:.3f}'.format(m_x) + "_lsig" + '{0:.2f}'.format(np.log10(sigma_p)) + ".txt"
headertxt = "mx [GeV]: " + str(m_x) + "\nsigma_p [cm^2]: " + str(sigma_p) + "\ndepth [m]: " + str(depth) + "\nloc: " + target
headertxt += "\nColumns: gamma/pi, v [km/s], f(v, gamma) [s/km]"

#Try different locations to save files, depending on where
#the code is being executed from
try:
    np.savetxt(fname, np.transpose([gamma_rep, vgrid.flatten(), fgrid.flatten()]), header=headertxt)
except:
    np.savetxt("../" + fname, np.transpose([gamma_rep, vgrid.flatten(), fgrid.flatten()]), header=headertxt)
    
