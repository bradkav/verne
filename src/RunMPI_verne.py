#!/usr/bin/env python
from mpi4py import MPI
from subprocess import call
import numpy as np

import sys

import argparse

"""
This script will create 'nprocs' parallel processes and on each one
calculate the speed distribution at the detector for a range of
values of gamma (using CalcVelDist.py).

It will then calculate the number of events at the detector using
CalcRate_CDMS.py or CalcRate_MPI.py

Run with:

    mpirun -np 16 python RunMPI_verne.py

Flags:
    -m_x: DM mass in GeV
    -target: "MPI" or "SUF" depending on location of detector
    -lsigstart, -lsigend: calculates for different values
                of log10(sigma) on a grid from lsigstart to lsigend
"""

#Parse the arguments!
parser = argparse.ArgumentParser(description='...')
parser.add_argument('-m_x','--m_x', help='DM mass in GeV', type=float,default = 1e5)
#parser.add_argument('-sigma_p','--sigma_p', help='DM-nucleon cross section, sigma_p in cm^2', type=float, required=True)
parser.add_argument('-target','--target', help='Location of detector: "MPI" or "SUF"', type=str, default="full")
parser.add_argument('-lsigstart', '--lsigstart', type=float, required=True)
parser.add_argument('-lsigend', '--lsigend', type=float,required=True) 
args = parser.parse_args()


comm = MPI.COMM_WORLD
#Get total number of MPI processes
nprocs = comm.Get_size()
#Get rank of current process
rank = comm.Get_rank()

#Sample between log10(sigma_p/cm^2) = -30.6...-27.6
#sig_list = np.logspace(-25.40, -22.40, nprocs)*(args.m_x/1e5)
sig_list = np.logspace(args.lsigstart, args.lsigend, nprocs)
sig = sig_list[rank]
#sig = 10**-27.9

#Output file labelled by r_np and exposure index
outfile = "outputs/out_"+args.target+"_lmx" + '{0:.1f}'.format(np.log10(args.m_x)) + "_lsig" + '{0:.2f}'.format(np.log10(sig))+".txt"

#Directory where the calc files are located
myDir = "/home/kavanagh/verne/"
cmd = "cd "+myDir+" ; python CalcVelDist.py "
cmd += "-m_x " + str(args.m_x)
cmd += " -sigma_p " + str(sig)
cmd += " -loc " + str(args.target)
cmd += " >> " + outfile

if (args.target == "SUF"):
    cmd += " ; python CalcRate_CDMS.py -m_x " + str(args.m_x)
elif (args.target == "MPI"):
    cmd += " ; python CalcRate_nucleus.py -m_x " + str(args.m_x)
    
cmd += " -sigma_p " + str(sig)
cmd += " >> " + outfile

sts = call(cmd,shell=True)
comm.Barrier()
