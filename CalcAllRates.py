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
from os import listdir
from subprocess import call


mstr = "5.0"
mx = 10**(float(mstr))

rootpath = "results/veldists/"

full_list = listdir(rootpath)
newlist = []

for f in full_list:
    if f.startswith("f_MPI_lmx" + mstr):
        newlist.append(f)

print " Number of files found:", len(newlist)
for f in newlist:
    inds = f.rfind("lsig")
    sigstr = f[inds+4:inds+10]
    sig = 10**(float(sigstr))
    cmd = "python CalcRate_nucleus.py -m_x " + str(mx) + " -sigma_p " + str(sig)
    call(cmd, shell=True)
    #print sigstr


#sts = call(cmd,shell=True)
