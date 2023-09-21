#!/bin/sh
#PBS -N pia2
#PBS -o vystup.txt
#PBS -l select=1:ncpus=4:mpiprocs=4
#PBS -M Frantisek.Stloukal@fs.cvut.cz
#PBS -m ae

cd $PBS_O_WORKDIR

module load GCC/8.3.0
module load OpenMPI/3.1.4-GCC-8.3.0

mpirun ./main
