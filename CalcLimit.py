import numpy as np
import WIMpy.DMUtils as DMU
from numpy import pi
from scipy.integrate import quad, dblquad, trapz, simps
import verne
from LabFuncs import *
from matplotlib.ticker import MultipleLocator, FormatStrFormatter
import utils
from scipy.special import erf
from timeit import default_timer as timer
from scipy.interpolate import interp2d, interp1d
from matplotlib import cm
from tqdm import tqdm
import argparse
import os.path

from matplotlib import pyplot as pl

exp = "nucleus"

if (exp == "CDMS"):
    loc = "SUF"
    N_lim = 34.96 #Upper limit for 27 observed events
    lmlist = ["1.0", "1.5", "2.0", "3.0", "4.0", "5.0"]
if (exp == "nucleus"):
    loc = "MPI"
    N_lim = 541.204 #Upper limit for 511 observed events
    lmlist = ["0.0", "1.0", "1.5", "2.0", "3.0", "4.0", "5.0"]


masses = np.zeros(len(lmlist))
limits = np.zeros(len(lmlist))

for i, lmstring in enumerate(lmlist):
    lsig, Ne = np.loadtxt("results/Nevents/N_" + loc+ "_lmx" + lmstring + ".txt", unpack=True)
    Ne = Ne[lsig.argsort()]
    lsig = np.sort(lsig) 
    
    Ne += (-lsig/1e20)

    #Generate an interpolating function
    lsig_interp = interp1d(np.log10(Ne), lsig)

    lsig_limit = lsig_interp(np.log10(N_lim))

    print lmstring, lsig_limit

    masses[i] = 10**float(lmstring)
    limits[i] = 10**lsig_limit

    pl.figure()
    pl.semilogy(lsig, Ne, "+")
    pl.semilogy(lsig, Ne, "-")
    pl.axvline(lsig_limit, linestyle='--', color='k')
    pl.show()
    
#Add on the final datapoint for very large masses
if (exp == "CDMS"):
    m_max = 1e15
if (exp == "nucleus"):
    m_max = 1e8

limits = np.append(limits, limits[-1]*m_max/masses[-1])
masses = np.append(masses, m_max)


    
np.savetxt("results/constraints/" + exp +"_thiswork.txt", zip(masses, limits))

#data = data[data[0,:].argsort()]
#data = np.sort(data, axis=1)

