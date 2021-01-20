#!/usr/bin/env python
from mpi4py import MPI
from subprocess import call
import numpy as np

import sys

import argparse


#Parse the arguments!
parser = argparse.ArgumentParser(description='...')
#parser.add_argument('-m_x','--m_x', help='DM mass in GeV', type=float,default = 1e5)
#parser.add_argument('-sigma_p','--sigma_p', help='DM-nucleon cross section, sigma_p in cm^2', type=float, required=True)
parser.add_argument('-target','--target', help='Target to propagate through. `earth`, `atmos`, `shield` or `full`.', type=str, default="full")
#parser.add_argument('-lsigstart', '--lsigstart', type=float, required=True)
#parser.add_argument('-lsigend', '--lsigend', type=float,required=True) 
parser.add_argument('-index', type=int)
args = parser.parse_args()

block_size = 10

comm = MPI.COMM_WORLD
#Get total number of MPI processes
nprocs = comm.Get_size()
#Get rank of current process
rank = comm.Get_rank()

mx_vals, sig_vals = np.loadtxt("data/param_grid.txt", unpack=True)


#Sample between log10(sigma_p/cm^2) = -30.6...-27.6
#sig_list = np.logspace(-25.40, -22.40, nprocs)*(args.m_x/1e5)
#sig_list = np.logspace(args.lsigstart, args.lsigend, nprocs)
#sig = sig_list[rank]
#sig = 10**-27.9

N_runs = nprocs*block_size

full_inds = np.arange(N_runs) + args.index*N_runs
this_inds = np.arange(block_size) + rank*block_size

index_list = full_inds[this_inds]

#current_m_vals = mx_vals[full_inds[this_inds]]
#current_sig_vals = sig_vals[full_inds[this_inds]]



#Directory where the calc files are located
myDir = "/home/kavanagh/verne/"
cmd = "cd "+myDir+" ; "

#for mx, sig in zip(current_m_vals, current_sig_vals):
for i in index_list:
    if (i < len(mx_vals)):
        mx = mx_vals[i]
        sig = sig_vals[i]
    
        #Output file labelled by r_np and exposure index
        outfile = "outputs/out_"+args.target+"_lmx" + '{0:.2f}'.format(np.log10(mx)) + "_lsig" + '{0:.2f}'.format(np.log10(sig))+".txt"
        
        cmd += "python3 src/CalcVelDist.py "
        cmd += " -m_x " + str(mx/1e3)
        cmd += " -sigma_p " + str(sig)
        cmd += " -loc " + str(args.target)
        cmd += " >> " + outfile + "; "

#calculate_rate = False
#
#if (calculate_rate):
#
#    if (args.target == "SUF"):
#        cmd += " ; python CalcRate_CDMS.py -m_x " + str(args.m_x)
#    elif (args.target == "MPI"):
#        cmd += " ; python CalcRate_nucleus.py -m_x " + str(args.m_x)
#    elif (args.target == "EDE"):
#        cmd += " ; python CalcRate_EDE.py -m_x " + str(args.m_x)
#    
#    cmd += " -sigma_p " + str(sig)
#    cmd += " >> " + outfile

print cmd
    
sts = call(cmd,shell=True)
comm.Barrier()
