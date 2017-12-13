#PBS -lnodes=1:cores16:ppn=16
#PBS -lwalltime=05:00:00

cd $HOME/verne/src

module load openmpi/gnu
#module load python/2.7.9

time mpiexec -np 16 python RunMPI_verne.py -m_x 1e3 -target SUF -lsigstart -29.77 -lsigend -29.62