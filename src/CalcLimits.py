import numpy as np
from scipy.interpolate import interp1d

import argparse
#Parse the arguments!
parser = argparse.ArgumentParser(description='...')
parser.add_argument('-exp','--exp', help='Calculate limits for which experiment? "CDMS" or "nucleus"', type=str, required=True)
args = parser.parse_args()
exp = args.exp


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

#Loop over the masses, calculating the upper limit
#From the Nevent files
for i, lmstring in enumerate(lmlist):
    lsig, Ne = np.loadtxt("../results/Nevents/N_" + loc+ "_lmx" + lmstring + ".txt", unpack=True)
    Ne = Ne[lsig.argsort()]
    lsig = np.sort(lsig) 
    
    Ne += (-lsig/1e20)

    #Generate an interpolating function
    lsig_interp = interp1d(np.log10(Ne), lsig)

    lsig_limit = lsig_interp(np.log10(N_lim))

    print lmstring, lsig_limit

    masses[i] = 10**float(lmstring)
    limits[i] = 10**lsig_limit
    
#Add on the final datapoint for very large masses
if (exp == "CDMS"):
    m_max = 1e15
if (exp == "nucleus"):
    m_max = 1e8

limits = np.append(limits, limits[-1]*m_max/masses[-1])
masses = np.append(masses, m_max)
    
np.savetxt("../results/constraints/" + exp +"_thiswork.txt", zip(masses, limits))


