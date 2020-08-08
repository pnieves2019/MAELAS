#!/bin/bash
#PBS -A OPEN-X-X
#PBS -q qprod
#PBS -l select=1:ncpus=24:mpiprocs=24:ompthreads=1
#PBS -l walltime=48:00:00
#PBS -N job_AA_BB
#PBS -j oe
#PBS -S /bin/bash

cd ${PBS_O_WORKDIR}
SCRDIR=/scratch/P_AA_BB

mkdir -p $SCRDIR
cd $SCRDIR || exit
cp -f -r $PBS_O_WORKDIR/* .
ml purge
ml VASP/5.4.4-intel-2017c-mkl=cluster
./vasp_0 >> log
exit
