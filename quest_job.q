#!/bin/bash
#MSUB -l nodes=2:ppn=20
##queues: short, normal, long
#MSUB -A b1027
#MSUB -q buyin
#MSUB -N job_name
#MSUB -m ae
#MSUB -M your_email@u.northwestern.edu
#MSUB -l walltime=hh:mm:ss

cd $PBS_O_WORKDIR
module load mpi/openmpi-1.6.3-intel2013.2

ulimit -s unlimited
nprocs=`wc -l $PBS_NODEFILE | awk '{ print $1 }'`
echo $PBS_NODEFILE

mpirun -n $nprocs /projects/b1027/VASPwannier90.5.3.5/vasp.5.3/vasp > output

