import numpy as np
from scipy.interpolate import interp1d

import MaxwellBoltzmann as MB
import argparse

import verne as verne
from matplotlib import pyplot as plt

try:
    from tqdm import tqdm
except:
    def tqdm(x):
        return x

def calcVelDist_full(m_x, sigma_p, loc, interaction, depth_in = 0):
    
    results_dir = "results/"


    if (interaction not in  ["SI", "SD", "hDP", "Millicharge"]):
        print(">Unknown interaction type <", interaction, ">...")
        exit()
        
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
    elif (loc == "full"):
        depth = depth_in
    
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
    
    v_th = 20e0 #Lowest speed to consider (don't go lower than v_th km/s, other the calculation of derivatives is messed up...)  
    
    def getVelDist(gamma):
        
        #It's useful to calculate the maximum final speed for particles reaching the detector from a given direction
        #in order not to waste time calculating the speed distribution for values of v_f which are impossible.
        print(">    Calculating maximum final speed...")
        a = 1.0
        b = 2*v_e*(-np.sin(gamma)*np.sin(np.pi-thetavals) + np.cos(gamma)*np.cos(np.pi-thetavals))
        c = v_e**2 - vesc**2
        v_initial_max = (-b + np.sqrt(b**2 - 4*a*c))/(2.0*a)
    
        #Calculate the maximum final speed as a function of incoming angle theta
        v_final_max = 0.0*v_initial_max
        for i in tqdm(range(Nvals)):
            v_final_max[i] = verne.calcVfinal_full(v_initial_max[i], thetavals[i],  depth, sigma_p, m_x, interaction, target)
    
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
        vlist = np.linspace(v_th, 0.9999999*vmax, 50)
        vlist = np.sort(vlist)
        f_final = 0.0*vlist
        for i in tqdm(range(len(vlist))):
            f_final[i] = verne.CalcF(vlist[i], gamma, depth, sigma_p, m_x, target, vfinal_interp, interaction=interaction)
    
        #Add on the final point
        vlist = np.append(vlist, vmax)
        f_final = np.append(f_final, 0.0)
    
        return vlist, f_final
        
        
    #Loop over gamma values
    N_gamma = 25
    Nv1 = 20 #Near the velocity threshold
    Nv2 = 78 #Everywhere else
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
    
