#!/opt/soft/python-2.7.3-1/bin/python
import argparse, os.path

# Author: Nick Wagner
# This file creates VASP job submission 
# files for the Carbon computing cluster 
# at Argonne National Lab. It requires 
# an INCAR with the static calculation 
# parameters named INCAR.static to exist

parser = argparse.ArgumentParser()
parser.add_argument("-n", action="store",
        dest="arguments", default="4_jobname_8",
        help='''Specify number of nodes, job name, and wall time
        in hours separated by underscores: e.g. Num_Name_Time. 
        4_jobname_8 is default''')
parser.add_argument("-jt", action='store', dest='jtype', default='static',
        help='''Job type. 'static', 'relax', and 'elast' are current options''')
parser.add_argument("-isif", action='store', dest='isif', default=3,
        help='''Relaxation scheme. 3 is default''')

email='nicholaswagner2018@u.northwestern.edu'

args = parser.parse_args()
isif = int(args.isif)
params = args.arguments.split("_")

# If POTCAR is not found, write one from POSCAR
if not os.path.isfile('POTCAR'):
    import wpot

# If KPOINTS.static is not found, write one
if not os.path.isfile('KPOINTS.static'):
    import kpointer

# If KPOINTS.bands is not found bu called for, write one
if not os.path.isfile('KPOINTS.bands') and jtype != 'elast':
    os.system("~/bin/kpath")

# Store ENCUT value for changing later
with open('INCAR.static', 'r') as inF:
    for line in inF:
        if 'ENCUT =' in line:
            encut = line.split()[2]

#PBS commands for header
with open('run_vasp.job', 'w') as f:
    f.write('#!/bin/bash\n')
    f.write('#PBS -N %s     # job name\n' % params[1])
    f.write('#PBS -o %s.out    #output file name\n' % params[1])
    f.write('#PBS -e %s.err    #error file name\n' % params[1])
    f.write('#PBS -l nodes=%s:ppn=8\n' % params[0])
    f.write('#PBS -l walltime=%s:00:00      #run time (hh:mm:ss)\n' % params[2])
    f.write('#PBS -A cnm*****  # account to charge\n')
    f.write('#PBS -M %s\n' % email)
    f.write('#PBS -m ea     # email when job finishes or aborts\n\n')
    f.write("cd $PBS_O_WORKDIR\n\n")

# Boilerplate for DOS and BS calculations
static_text = '''sed -i "s/NSW = .*/NSW = 0 # number of ionic steps/" INCAR.static
sed -i "s/LCHARG = .FALSE./LCHARG = .TRUE./" INCAR.static
sed -i "s/ISMEAR = .*/ISMEAR = -5/" INCAR.static
sed -i "s/LREAL = .*/LREAL = .FALSE./" INCAR*
cp INCAR.static INCAR
cp KPOINTS.static KPOINTS
mpirun -machinefile $PBS_NODEFILE -np $PBS_NP vasp > out.static
mv OUTCAR OUTCAR.static
mv vasprun.xml %s_dos.xml
        
sed "s/ICHARG = .*/ICHARG = 11/" INCAR.static > INCAR.bands
sed -i "s/ISMEAR = .*/ISMEAR = 0/" INCAR.bands
sed -i "s/.*SIGMA = .*/SIGMA = 0.1/" INCAR.bands
sed -i "s/.*LCHARG = .*/LCHARG = .FALSE./" INCAR.bands
cp KPOINTS.bands KPOINTS
cp INCAR.bands INCAR
mpirun -machinefile $PBS_NODEFILE -np $PBS_NP vasp > out.bands
mv vasprun.xml %s_BS.xml
mv OUTCAR OUTCAR.bands
rm WAVECAR\n'''


# Boilerplate for ISIF relaxations
relax_text = '''cp KPOINTS.static KPOINTS
cp POSCAR POSCAR.orig
sed "s/ISMEAR = .*/ISMEAR = 0/" INCAR.static > INCAR.is%s.ib1
sed -i "s/.*SIGMA = .*/SIGMA = 0.1/" INCAR.is%s.ib1
sed -i "s/.*LCHARG = .*/LCHARG = .FALSE./" INCAR.is%s.ib1
sed -i "s/NSW = .*/NSW = 40/" INCAR.is%s.ib1
sed -i "s/LREAL = .*/LREAL = .FALSE./" INCAR*
sed -i "s/.*ISIF = .*/ISIF = %s/" INCAR.is%s.ib1
sed -i "s/.*IBRION = .*/IBRION = 1/" INCAR.is%s.ib1
sed -i "s/.*POTIM = .*/POTIM = 0.2/" INCAR.is%s.ib1
sed -i "s/.*EDIFFG = .*/EDIFFG = -0.0005/" INCAR.is%s.ib1
sed -i "s/ENCUT = .*/ENCUT = %s/" INCAR.is%s.ib1
sed "s/IBRION = 1/IBRION = 2/" INCAR.is%s.ib1 > INCAR.is%s.ib2
for it in 1 2 3
do
cp INCAR.is%s.ib2 INCAR
mpirun -machinefile $PBS_NODEFILE -np $PBS_NP vasp > out.is%s.ib2.$it
mv OUTCAR OUTCAR.is%s.ib2.$it
cp CONTCAR CONTCAR.is%s.ib2.$it
cp  CONTCAR.is%s.ib2.$it POSCAR
cp INCAR.is%s.ib1 INCAR
mpirun -machinefile $PBS_NODEFILE -np $PBS_NP vasp > out.is%s.ib1.$it
mv OUTCAR OUTCAR.is%s.ib1.$it
cp CONTCAR CONTCAR.is%s.ib1.$it
cp CONTCAR POSCAR
done\n\n'''

# Elast text is similar to static but we don't calculate band structure or save the CHGCAR
# We also leave the vasprun alone rather than renaming it.
elast_text = '''sed -i "s/NSW = .*/NSW = 0 # number of ionic steps/" INCAR.static
sed -i "s/LCHARG = .FALSE./LCHARG = .TRUE./" INCAR.static
sed -i "s/ISMEAR = .*/ISMEAR = -5/" INCAR.static
sed -i "s/LREAL = .*/LREAL = .FALSE./" INCAR*
cp INCAR.static INCAR
cp KPOINTS.static KPOINTS
mpirun -machinefile $PBS_NODEFILE -np $PBS_NP vasp > out.static
mv OUTCAR OUTCAR.static
rm CHGCAR WAVECAR\n'''

#Commands for calculating DOS and band structure
if args.jtype == 'static':
    with open('run_vasp.job', 'a') as f:
        f.write(static_text % (params[1], params[1]))


# Commands for ionic relaxation
if args.jtype == 'relax':
    rel_encut = str(int(encut)*0.9)
    with open('run_vasp.job', 'a') as f:
        # Format relax_text with isif number and encut number
        a = [isif]*9
        a.append(rel_encut)
        a.extend([isif]*12)
        f.write(relax_text % tuple(a))
        f.write(static_text % (params[1], params[1]))
        
if args.jtype == 'elast'
    rel_encut = str(int(encut)*0.9)
    isif = 2
    with open('run_vasp.job', 'a') as f:
        # Format relax_text with isif number and encut number
        a = [isif]*9
        a.append(rel_encut)
        a.extend([isif]*12)
        f.write(relax_text % tuple(a))
        f.write(elast_text % (params[1], params[1]))
