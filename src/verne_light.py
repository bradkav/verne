"""
verne_light.py - Code for calculating the Earth stopping effect, primarily for light (MeV) Dark Matter.


Last updated: 24/10/2024
Contact: Bradley Kavanagh, bradkav@gmail.com

"""

import numpy as np
from scipy.integrate import simps

import scipy.special
import os
import sys

import MaxwellBoltzmann as MB
import density
from density import R_E, h_A

try:
    from tqdm import tqdm
except:
    tqdm = lambda x: x


#--------------------
#Theta = 0 is directly from BELOW, angles in radians
#Gamma = 0 is mean DM flux directly from BELOW, angles in radians

#Densities are in number/cm^3
#Distances in m
#--------------------

m_e = 511e-3 #Electron mass in GeV
alpha = 1/137.03 #Fine structure constant
q_screen = alpha*m_e #Electron screening momentum in GeV

#--------------------
#Integration parameters
N_THETA = 201

#--------------------

global X_grid
X_grid = None

#Indices of different elements in the Earth and atmosphere
isovals = range(density.Niso_full)

#Grid of theta values for interpolation
thetas_full = np.linspace(0, np.pi, 10001)

    
#Path length [m], as measured from the top of the atmosphere to the detector 
#(at 'depth' m underground)
def pathLength(depth, theta):
    r_det = R_E - depth
    return +np.cos(theta)*r_det + np.sqrt((-np.cos(theta)*r_det)**2 - (r_det**2 - (R_E+h_A)**2))
   
#Path length [m], as measured from the Earth's surface to the detector 
#(at 'depth' m underground)
def pathLength_Earth(depth, theta):
    r_det = R_E - depth
    return +np.cos(theta)*r_det + np.sqrt((-np.cos(theta)*r_det)**2 - (r_det**2 - (R_E)**2))
    
#Maximum recoil energy (in keV)
def ERmax(mX, mA, v):
    mu = mX*mA*1.0/(mX+mA)
    return (1e6/(3e5*3e5))*2*(mu*v)**2/mA
    
#Generate a grid of theta values, over the relevant range for a given 
#velocity distribution and isodetection angle
def generate_theta_grid(v, gamma):
   fudge = 1e-6
    
   A = (MB.vesc**2 - v**2 - MB.ve**2)/(2*(v + fudge)*MB.ve)
   A = np.clip(A, -1, 1)

   #Range of thetas set by requiring |v - v_e|^2 < v_esc^2
   theta_min = (1 + fudge)*np.clip(np.arccos(A) + gamma - np.pi, 0, np.pi)
   if theta_min == 0:
       theta_min = fudge
   theta_max = (1 - fudge)*np.clip(np.pi - (np.arccos(A) - gamma), 0, np.pi)
   
   x_list = np.linspace(0, np.pi, N_THETA)
   thetas = theta_min + 0.5*(np.cos(x_list[::-1]) + 1)*(theta_max - theta_min)
   
   return thetas
    

def Tabulate_Column_Density(depth, target="full"):
    global X_grid
    
    X_grid = np.zeros((len(isovals), len(thetas_full))) 

    for i, theta in enumerate(tqdm(thetas_full, desc="Tabulating Column Density")):
        if (target == "atmos"):
            d1 = 0
            d2 = pathLength(depth, theta) - pathLength_Earth(depth, theta)
            d_vec = np.linspace(d1, d2, 100)
        elif (target == "earth"):
            d1 = pathLength(depth, theta) - pathLength_Earth(depth, theta)
            d2 = pathLength(depth, theta)
            d_vec = np.linespace(d1, d2, 100)
        elif (target == "full"):
            d1 = 0
            ds = pathLength(depth, theta) - pathLength_Earth(depth, theta)
            d2 = pathLength(depth, theta)
            d_vec = np.append(np.linspace(d1, ds, 20), np.linspace(ds*1.01, d2, 80))
    
        r_vec = radius(d_vec, theta, depth)
        
        for j in range(len(isovals)):
            n_j = density.dens_interp(isovals[j], r_vec)
            X_grid[j,i] = simps(n_j, d_vec)
            
            
def Calc_p_trans_grid(thetas, v, sigma_p, m_x, target, interaction="SI"):
    p_trans = 0.0*thetas

    for i in range(len(thetas)):
        params = [sigma_p, m_x, interaction, target, v]
   
        #Distance to the detector divided by mean free path (averaged over trajectory)
        L1_over_lmbda_avg = inv_mean_free_path_avg_r(thetas[i], params)
        
        ### Transmisson probability
        expL1 = np.exp(-L1_over_lmbda_avg)
        P2_L1 = L1_over_lmbda_avg/2*expL1 - expL1/4 + expL1**3/4
        P0_L1 = expL1
        
        p_trans[i] = (P0_L1 + P2_L1)
    
    return p_trans


#Calculate the final speed distribution at the detector
def CalcF_transmitted(v, gamma, sigma_p, m_x, target, interaction="SI"):
    
    if (X_grid is None):
        raise Exception("Column density has not been tabulated. Call `verne.Tabulate_Column_Density(depth, target)` before calculating velocity distributions.")
    
    thetas = generate_theta_grid(v, gamma)
    
    p_trans = Calc_p_trans_grid(thetas, v, sigma_p, m_x, target, interaction)
    
         
    #f_full = 0.0*thetas
    #for i, theta in enumerate(thetas):
    #    f_full[i] = f_integrand_full(v, theta, gamma)
    f_full = f_integrand_full(v, thetas, gamma)
    f_full = np.nan_to_num(f_full)
    
    fint_2order = f_full*p_trans

    #Integrate with Simpson's rule
    return simps(fint_2order, thetas)
 

def Calc_p_refl_grid(thetas, v, sigma_p, m_x, target, interaction="SI"):

    p_refl = 0.0*thetas 
    
    for i in range(len(thetas)):
        params = [sigma_p, m_x, interaction, target, v]
        
        #Particles reflecting passed the detector
        #Distance _after passing_ the detector divided by mean free path (averaged over trajectory)
        L2_over_lmbda_avg = inv_mean_free_path_avg_r(np.pi - thetas[i], params)
        expL2 = np.exp(-L2_over_lmbda_avg)
        
        ### Reflection probability
        P2_L2 = L2_over_lmbda_avg/2*expL2 - expL2/4 + expL2**3/4
        P0_L2 = np.exp(-L2_over_lmbda_avg)
        
        #Path before the detector
        #Distance to the detector divided by mean free path (averaged over trajectory)
        L1_over_lmbda_avg = inv_mean_free_path_avg_r(thetas[i], params)
        expL1 = np.exp(-L1_over_lmbda_avg)
        P2_L1 = L1_over_lmbda_avg/2*expL1 - expL1/4 + expL1**3/4
        P0_L1 = expL1

        p_refl[i] = (1-P0_L2-P2_L2)*(P0_L1+P2_L1)
    
    return p_refl

#Calculate the reflected component of the distribution at the detector
def CalcF_reflected(v, gamma, sigma_p, m_x, target, interaction="SI"):
    
    if (X_grid is None):
        raise Exception("Column density has not been tabulated. Call `Calc_Column_Density(depth, target)` before calculating velocity distributions.")
    
    thetas = generate_theta_grid(v, gamma)
    p_refl = Calc_p_refl_grid(thetas, v, sigma_p, m_x, target, interaction)
             
    #f_full = 0.0*thetas
    #for i, theta in enumerate(thetas):
    #    f_full[i] = f_integrand_full(v, theta, gamma)

    f_full = f_integrand_full(v, thetas, gamma)

    fint_2order = f_full*p_refl

    #Integrate with Simpson's rule
    return simps(fint_2order, thetas)


#Integrand for calculating the final speed distribution at the detector
def f_integrand_full(v, theta, gamma):
    return np.sin(theta)*(v**2)*MB.calcf_integ(v, theta, gamma)
 

#Calculate the distance of a point from the centre of the Earth
#The point is defined by:
#   - theta, the angle of the trajectory
#   - depth,the detector depth
#   - D, the distance along the trajectory, starting at the top of the atmosphere
def radius(D, theta, depth):
    r_det = R_E - depth 
    return np.sqrt((R_E+h_A)**2 + D**2 + 2*D*(r_det*np.cos(theta) - pathLength(depth, theta)))

#Correction of the cross-section due to charge screening
def sigma_correction(Emax, Z, A, interaction):
    m_A = A*931.5*1e3 # keV
    a0 = 52917.7249 #fm
    q1 = np.sqrt(2*m_A*Emax) # keV
    keVfm_1 = 1e-12/1.97e-7 #1/keV/fm
    q = q1*keVfm_1 #fm^-1
    a = 1/4*(9*np.pi**2/(2*Z))**(1/3)*a0#/keVfm_1
    aq = a*q
    
    if interaction.lower() == "hm":
        p_refl = 7/8
        return p_refl*( 1 + 1/(1+aq**2) - 2/aq**2*np.log(1+aq**2) )
        
    elif interaction.lower() == "ulm":
        me = 0.511e3 #keV
        qref = me*1/137*keVfm_1 #me * alpha #fm^-1
        p_refl = 1/2
        return p_refl*a**4*qref**4/(1+aq**2)


### Calculates the inverse of the mean free path for an array of n
def inv_mean_free_path_vec(sigma_p, m_X, n, A, mediator=False,interaction="SI"):
    m_p = 0.9315 #Proton mass
    if mediator:
        m_A = m_p*mediator[1]
    else:
        m_A = m_p*A
    mu_A = m_A*m_X/(m_A + m_X)
    mu_p = m_p*m_X/(m_p + m_X)
    if mediator:
        return ( n*A**2*mu_A**2 )*sigma_p*sigma_correction(ERmax(m_X,0.9315*mediator[1], mediator[0]), A, mediator[1],interaction)/(1e-2*mu_p**2) 
    else:
        #For very light DM, reflection probability with scattering is 50%
        p_refl = 1/2
        return p_refl*(n*sigma_p*A**2*mu_A**2 )/(1e-2*mu_p**2) 


### Calculates the averaged inverse of the mean free path for an array of n
### over all the isotopes
def inv_mean_free_path_avg_r(theta, params):

    sigma_p, m_x, interaction, target, v = params

    if (target == "atmos"):
        isovals = [8,9]
    elif (target == "earth"):
        isovals = range(density.Niso)
    else:
        isovals = range(density.Niso_full)
    
    #r = radius(D, theta, depth)
    if (interaction.lower() == "SI".lower()):
        C = density.Avals
    elif (interaction.lower() == "hm" or interaction.lower() =="ulm"):
        C = density.Zvals 
 
    #Loop over the relevant isotopes
    inv_mean_free_path_avg = 0.0
    for i in isovals:
        X = np.interp(theta, thetas_full, X_grid[i,:])
        #n_i = density.dens_interp(i, r)
        #mean_free_path_avg += mean_free_path_vec(sigma_p, m_x, n_i, Avals[i])

        if (interaction.lower() == "hm" or interaction.lower() =="ulm"):
            inv_mean_free_path_avg += inv_mean_free_path_vec(sigma_p, m_x, X, C[i], mediator=[v, density.Avals[i]],interaction=interaction)
        else:
            inv_mean_free_path_avg += inv_mean_free_path_vec(sigma_p, m_x, X, C[i])
    return inv_mean_free_path_avg


