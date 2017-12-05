#PBS -lnodes=1:cores16:ppn=16
#PBS -lwalltime=03:00:00

cd $HOME/verne/

module load openmpi/gnu
#module load python/2.7.9

time mpiexec -np 16 python RunMPI_verne.py -m_x 1e5 -target full