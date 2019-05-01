# Job-Scripts
Scripts for submitting VASP calculations to Quest, Carbon, and Stampede supercomputers.

Check the VASP user guide for info on flags. 
https://cms.mpi.univie.ac.at/vasp/vasp/vasp.html

For those running jobs on SLURM schedulers, here are some aliases that may be useful:
* `alias sq='squeue -l'` #look at all jobs
* `alias mtdq='squeue -A <alloc code>'` #look at jobs of a specific allocation
* `alias myq='squeue -u <username>'` #look at your jobs
* `alias st='squeue --start -u <username>'` #get estimated start times for your pending jobs
