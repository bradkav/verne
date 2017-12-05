#!/usr/bin/env python
from mpi4py import MPI
from subprocess import call
import sys

import argparse

"""
This script will create 'nprocs' parallel instances of 
CalcDisc-vs-Exposure.py and calculate the significance
for a different exposure on each instance.

For a given mass and ensemble, there are 32 exposures to
calculate for, so this code is meant to be run in batches. 
The number of the batch is specified in arr_index (starting
from arr_index=1). The code will then automatically calculate 
the appropriate index which labels the exposure to be used,
based on the batch number and the number of processes in 
the current batch.

See CalcDisc-vs-Exposure.py for information on what the other
command line parameters should be.

To execute with 16 MPI processes, you would have something
like the following in your submission script:

mpirun -np 16 python RunMPI_exposure.py
"""

#Parse the arguments!
parser = argparse.ArgumentParser(description='...')
parser.add_argument('-m_x','--m_x', help='DM mass in GeV', type=float,default = 1e5)
#parser.add_argument('-sigma_p','--sigma_p', help='DM-nucleon cross section, sigma_p in cm^2', type=float, required=True)
parser.add_argument('-target','--target', help='Target to propagate through. `earth`, `atmos`, `shield` or `full`.', type=str, default="full")
args = parser.parse_args()


comm = MPI.COMM_WORLD
#Get total number of MPI processes
nprocs = comm.Get_size()
#Get rank of current process
rank = comm.Get_rank()

#Sample between log10(sigma_p/cm^2) = -30.6...-27.6
lsig_list = np.logspace(-30.6, -27.6, nprocs)
lsig = lsig_list[rank]

#Output file labelled by r_np and exposure index
outfile = "outfiles/f_lmx" + '{0:.1f}'.format(np.log10(m_x)) + ".txt"

#Directory where the calc files are located
myDir = "/home/kavanagh/verne/"
cmd = "cd "+myDir+" ; python CalcVelDist.py "
cmd += "-m_x " + str(args.m_x)
cmd += "-sigma_p " + str(lsig)
cmd += "-target " + str(args.target)
cmd += ">> " + outfile

sts = call(cmd,shell=True)
comm.Barrier()
