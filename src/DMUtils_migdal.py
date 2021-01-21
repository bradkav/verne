import numpy as np

from scipy.interpolate import interp1d
from scipy.integrate import quad, simps
from scipy.special import erf

import os

#Some parameters
rho_x = 0.3 #GeV cm^-3

file_dir = os.path.dirname(os.path.realpath(__file__))

def load_probabilities(target, n, l):
    
    readdata = False
    readQN = False
    
    Evals = []
    pvals = []

    
    filename = file_dir + "/../data/" + target + ".dat"
    with open(filename) as f:
        for line in f:
            linestr = line.split()
            #print(linestr[0])
            if (readdata == True):
                if ('Principal' in linestr):
                    break
                    #return np.array(Evals), np.array(pvals)
                    
                if ('Electron' not in linestr):
                    Evals.append(float(linestr[0]))
                    pvals.append(float(linestr[1]))
            
            if (readQN):
                n_current = int(linestr[0])
                l_current = int(linestr[1])
                #print(n_current, l_current)
                readQN = False
                if (n_current == n and l_current == l):
                    readdata = True
                    #print(n_current, l_current)
        
            if ('Principal' in linestr):
                #print("Here")
                readQN = True
    
    if (readdata == False):
        Evals = np.logspace(np.log10(1), np.log10(7e4),251)
        pvals = np.zeros(251)
    #eV, eV^-1
    return np.array(Evals), np.array(pvals)

def calc_pinterp(target, n, l):
    
    Evals, pvals = load_probabilities(target, n, l)
    return interp1d(Evals*1e-3, pvals*1e3, bounds_error=False, fill_value=0.0)

p_trans = [ [calc_pinterp("Ge", n, l) for l in [0,1,2]] for n in [1,2,3,4]]


def calc_ptrans(n,l, qe, E):
    return qe**2*p_trans[n-1][l](E)

def calc_qesq(E_R, A):
    m_A = 0.9315e9*A
    m_e = 511e3
    return 2*m_e**2*E_R*1e3/m_A

def calc_ptrans_full(n,l, E_R, E_e, A):
    qesq = calc_qesq(E_R, A)
    return qesq*p_trans[n-1][l](E_e)

#In km/s
def vmin_migdal(E_R, m_x, A, deltaE):
    m_A = 0.9315*A
    mu = m_x*m_A/(m_x + m_A)
    return 3e5*(m_A*E_R*1e-6 + mu*deltaE*1e-6)/(mu*np.sqrt(2*m_A*E_R*1e-6))

#-----------------------------------------------------------
# Standard Helm Form Factor for SI scattering
def calcSIFormFactor(E, m_N, old=False):

    #Define conversion factor from amu-->keV
    amu = 931.5*1e3

    #Convert recoil energy to momentum transfer q in keV
    q1 = np.sqrt(2*m_N*amu*E)

    #Convert q into fm^-1
    q2 = q1*(1e-12/1.97e-7)
    
    #Calculate nuclear parameters
    s = 0.9
    a = 0.52
    c = 1.23*(m_N**(1.0/3.0)) - 0.60
    R1 = np.sqrt(c*c + 7*np.pi*np.pi*a*a/3.0 - 5*s*s)
    
    if (old):
        R1 = np.sqrt((1.2**2)*m_N**(2.0/3.0) - 5)
 
    x = q2*R1
    J1 = np.sin(x)/x**2 - np.cos(x)/x
    F = 3*J1/x
    
    formfactor = (F**2)*(np.exp(-(q2*s)**2))
    #formfactor[E < 1e-3] = 1.0
    return formfactor

#-----------------------------------------------------------
# Reduced mass - input A as nucleon number and m_x in GeV
def reduced_m(A, m_x):
    m_A = 0.9315*A
    return (m_A * m_x)/(m_A + m_x)


def E_nl(target, n ,l):
    if (target == "Ge"):
        if (n == 1 and l == 0):
            return 1.1e4
        if (n == 2 and l == 0):
            return 1.4e3
        if (n == 2 and l == 1):
            return 1.2e3
        if (n == 3 and l == 0):
            return 1.7e2
        if (n == 3 and l == 1):
            return 1.2e2
        if (n == 3 and l == 2):
            return 3.5e1
        if (n == 4 and l == 0):
            return 1.5e1
        if (n == 4 and l == 1):
            return 6.5

conversion_factor = 1e-6*(1.79e-27)**-1*(3e5)**2*1e5*(60*60*24)
#Events for keV per kg per day per km/s
def d2RdEdv(E_R, v, m_x, sigma_p, A, speeddist=None, SD=None):
    if (SD == None):
        int_factor = sigma_p*A**2*calcSIFormFactor(E_R, A)
    else:
        J_Ge = 9/2
        if (SD == "p"):
            spin_sq = S_p**2
        elif (SD == "n"):
            spin_sq = S_n**2
        #print(SD)
        int_factor = (4/3)*((J_Ge+1)/J_Ge)*spin_sq*FormFactor_SD(E_R, SD)*sigma_p
    mu = reduced_m(1.0, m_x)
    prefactor = 0.5*rho_x/(m_x*mu**2)
    
    if (speeddist == None):
        speeddist = calcf_SHM
    
    #Does all the information about v factorise?
    #*NO*, because vmin depends on the ionisation process...
    return conversion_factor*prefactor*int_factor*speeddist(v)/v

nl_pair = [[1,0],
           [2,0],
           [2,1],
           [3,0],
           [3,1],
           [3,2],
           [4,0],
           [4,1]]

def d3R(E_R, E_EM, v, m_x, sigma_p, A, n_sep, speeddist=None, SD=None):
    dR = d2RdEdv(E_R, v, m_x, sigma_p, A, speeddist, SD)
    Z = 0
    for i in range(len(nl_pair)):
        n = nl_pair[i][0]
        if (n_sep == n or n_sep == -1):
            l = nl_pair[i][1]
            E_e = E_EM - E_nl("Ge", n,l)/1e3
            Z += calc_ptrans_full(n,l, E_R, E_e, A)
    return dR*Z/(2*np.pi)

def d2RdERdEe(E_R, E_EM, m_x, sigma_p, A, n_sep=-1, speeddist=None, SD=None, force_quad=False):
    integ = lambda v: d3R(E_R, E_EM, v, m_x, sigma_p, A, n_sep, speeddist, SD)
    vmin = vmin_migdal(E_R, m_x, A, E_EM)
    if (force_quad):
        res = quad(integ, vmin, 800)[0]
    else:
        vlist = np.linspace(vmin, 800, 50)
        integ_list = np.vectorize(integ)(vlist)
        res = simps(integ_list, vlist)
    return res

#------------
# Speed distribution

vesc = 533.0
sigmav = 156.0
ve = 232.0
# Nesc - normalisation constant
Nesc = (erf(vesc/(np.sqrt(2.0)*sigmav)) - np.sqrt(2.0/np.pi)*(vesc/sigmav)*np.exp(-vesc**2/(2.0*sigmav**2)))

def calcf_SHM(v):
    
    aplus = np.minimum((v+ve), v*0.0 + vesc)/(np.sqrt(2)*sigmav)
    aminus = np.minimum((v-ve), v*0.0 + vesc)/(np.sqrt(2)*sigmav)
    
    f = np.exp(-aminus**2) - np.exp(-aplus**2)


    return v*f/(np.sqrt(2*np.pi)*sigmav*ve*Nesc)




#-----------------
# Spin-dependent interactions
#-----------------
S_p = 0.031
S_n = 0.439

coeff_S00 = np.array([ 0.215608 , - 0.578786 , 0.698020, - 0.372000 , 0.107576 , - 0.0182408 , 0.00217108 , - 2.07981e-4, 1.65907e-5, - 5.95664e-7 ])
coeff_S11_min = np.array([  0.0743728 , - 0.233814 , 0.341725 , - 0.259024 , 0.121206 , - 0.0371226 , 0.00741080, - 9.02610e-4 , 5.81933e-5, - 1.38557e-6 ])
coeff_S11_max = np.array([  0.120045 , - 0.384157 , 0.559728 , - 0.415686 , 0.188412 , - 0.0568025, 0.0120204 , - 0.00175855 , 1.59975e-4 , - 6.66472e-6 ])
coeff_S01_min = np.array([  - 0.321836 , 0.950136 , - 1.27413 , 0.831035 , - 0.323769 ,  0.0831244 , - 0.0151542 , 0.00193259 , - 1.55025e-4, 5.68777e-6 ])
coeff_S01_max = np.array([  - 0.253289 , 0.739394 , - 0.993188 , 0.659953 , - 0.269522 , 0.0745897, - 0.0144162 , 0.00181542 , - 1.29365e-4, 3.77020e-6 ])

coeff_Sp_min = np.array([ 0.0138433 , - 0.0138982 , - 0.00961825 , 0.0275620 , - 0.0101577 , - 0.00235492 , 0.00246030, - 6.53041e-4 , 7.84526e-5, - 3.61078e-6 ])
coeff_Sp_max = np.array([ 0.0366954 , - 0.0733258 , 0.0471313 , 0.0281229 , - 0.0405538 , 0.0196085 , - 0.00515247 , 8.06626e-4, - 6.95571e-5, 2.63102e-6 ])
coeff_Sn_min = np.array([ 0.543270 , - 1.55198 , 2.03269 , - 1.28990 , 0.496419 , - 0.128347 , 0.0232676 , - 0.00274482 , 1.81026e-4 , - 4.56383e-6 ])
coeff_Sn_max = np.array([ 0.657509 , - 1.91400 , 2.53820 , - 1.63488 , 0.639763 , - 0.171656 ,  0.0345442 , - 0.00504185 , 4.64828e-4 , - 1.93402e-5 ])

b = 2.1058 #fm

def getSSF(u, coeff):
    return np.exp(-u)*np.polyval(coeff[::-1], u)

def SSF_p_min(u):
    return getSSF(u, coeff_Sp_min)

def SSF_p_max(u):
    return getSSF(u, coeff_Sp_max)

def SSF_p(u):
    return 0.5*(getSSF(u, coeff_Sp_min) + getSSF(u, coeff_Sp_max))

#----------------

def SSF_n_min(u):
    return getSSF(u, coeff_Sn_min)

def SSF_n_max(u):
    return getSSF(u, coeff_Sn_max)

def SSF_n(u):
    return 0.5*(getSSF(u, coeff_Sn_min) + getSSF(u, coeff_Sn_max))

#-----------------

def FormFactor_SD(E_R, nucleon="p"):
    m_N = 73.0
    
    #Define conversion factor from amu-->keV
    amu = 931.5*1e3

    #Convert recoil energy to momentum transfer q in keV
    q1 = np.sqrt(2*m_N*amu*E_R)
    
    #Convert q into fm^-1
    q2 = q1*(1e-12/1.97e-7)
    
    u = 0.5*(q2*b)**2
    
    if (nucleon == "p"):
        return SSF_p(u)/SSF_p(0)
    elif (nucleon == "n"):
        return SSF_n(u)/SSF_n(0)
    
#------------------