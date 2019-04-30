#!/bin/bash
#SBATCH -A b1027        # which account to debit hours from
#SBATCH -J Test               # job name
#SBATCH -o myjob.o%j           # output and error file name (%j expands to jobID) 
#SBATCH -e myjob.e%j           # output and error file name (%j expands to jobID) 
#SBATCH -N 2
#SBATCH -n 40                  # total number of mpi tasks requested
#SBATCH -p buyin              # queue (partition) -- normal, development, etc.
#SBATCH -t 2:00:00            # wall time (hh:mm:ss)
#SBATCH --mail-user=YourName@u.northwestern.edu 
#SBATCH --mail-type=begin      # email when job starts
#SBATCH --mail-type=end        # email when job ends
# Uses default modules: intel/2016.0, mpi/openmpi-1.6.3-intel2013.2  
mpirun -n $SLURM_NTASKS /projects/b1027/VASPmod.5.4.4/vasp_std > out
