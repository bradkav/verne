import numpy as np
import WIMpy.DMUtils as DMU
from numpy import pi
from scipy.integrate import quad


# Some random functions...

def q(A, E_R):
    mN = 1e3*A*0.9315 #in MeV
    
    return np.sqrt(2*mN*E_R*1e-3) #in MeV

def mu(mX, mA):
    return mX*mA*1.0/(mX+mA)

def ERmax(mX, mA, v):
    return 2*(mu(mX, mA)*v)**2/mA
    
#--------------------

A = 28

m_N = A*0.9315
m_x = 1e5

nA = 1e23 # in cm^-3
sigma_p = 1e-27 # in cm^2



#Neglecting any form factors for now!
sigma_A = (DMU.reduced_m(A, m_x)/DMU.reduced_m(1.0, m_x))**2*A**2*sigma_p
correction = quad(lambda x: DMU.calcSIFormFactor((1e6/(3e5*3e5))*x*ERmax(m_x, m_N, 220), A), 0, 1)[0]
print sigma_A, correction

Lambda = 1e-5*1.0/(nA*sigma_A) # in km



def step(x, v, dt):
    magv = np.sqrt(v[0]**2 + v[1]**2 + v[2]**2)
    dx = magv*dt
    
    if (x[2] > 100*2*6371.0):
        vp = v
        xp = x + v*dt
        return (xp, vp)
    
    #Check to see if a scatter happened
    if (np.exp(-dx/Lambda) < 0.9):
        print " Warning: dp > 0.1. Consider reducing dt..."
    if (np.random.rand(1) < (1 - np.exp(-dx/Lambda))):
        #print " Scatter!"
        #Assume SI scattering - still neglecting form factors!!!
        ER = ERmax(m_x, m_N, magv)*np.random.rand(1)
        Ei = 0.5*m_x*magv**2

        Ef = Ei - ER
        vf = np.sqrt(2*Ef/(m_x))
        #print (vf-magv)/magv
        q = np.sqrt(2*m_N*ER)
        ctp = (magv - q**2/(2*m_x*mu(m_x, m_N)*magv))/vf
        ctp = np.clip(ctp, -1, 1)
        tp = np.arccos(ctp)
        
        t1 = np.arccos(v[2]/magv)
        p1 = np.arctan2(v[1],v[0])
        
        pr = 2*pi*np.random.rand(1)
        vint = v*0.0
        vint[0] = vf*np.sin(tp)*np.cos(pr)
        vint[1] = vf*np.sin(tp)*np.sin(pr)
        vint[2] = vf*np.cos(tp)
        
        vp = 0.0*v
        
        #Rotate into correct coordinate frame
        vp[0] = np.cos(p1)*np.cos(t1)*vint[0] - np.sin(p1)*vint[1] +np.cos(p1)*np.sin(t1)*vint[2]
        vp[1] = np.sin(p1)*np.cos(t1)*vint[0] + np.cos(p1)*vint[1] +np.sin(p1)*np.sin(t1)*vint[2]
        vp[2] = -np.sin(t1)*vint[0] + np.cos(t1)*vint[2]
        
        xp = x + vp*dt
        #print Ef
        
    else:
        vp = v
        xp = x + v*dt
    
    return (xp, vp, dx)
    