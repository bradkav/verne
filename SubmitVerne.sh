#!/bin/bash
 
#SBATCH -N 1 --ntasks-per-node=16  
#SBATCH -t 06:00:00

#SBATCH -p normal

#SBATCH -o slurm_output/slurm-%j.out # STDOUT
#SBATCH -e slurm_output/slurm-%j.err # STDERR

cd $HOME/verne/

#module load openmpi/gnu
#module load python/2.7.9

module load pre2019

module unload GCCcore
module load Python/2.7.12-intel-2016b
module load slurm-tools

export SLURM_CPU_BIND=none

#time mpirun -np 16 python2.7 RunMPI_verne.py -m_x 0.200 -target EDE -lsigstart -25.50 -lsigend -24.00
time mpirun -np 16 python2.7 RunMPI_verne.py -target MOD -index $1
