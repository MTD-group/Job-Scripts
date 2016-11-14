#!/home/naw356/bin/python
import argparse, os.path

# If POTCAR is not found, write one from POSCAR
if not os.path.isfile('POTCAR'):
    import wpot

parser = argparse.ArgumentParser()
parser.add_argument("-n", action="store",
        dest="arguments", default="b1029_2_jobname_8",
        help='''Specify account number, number of nodes, job name, and wall time
        in hours separated by underscores: e.g. Acct_NodeNum_Name_Time''')
parser.add_argument("-jt", action='store', dest='jtype', default='static',
        help='''Job type. 'static' and 'relax' are current options''')
parser.add_argument("-isif", action='store', dest='isif', default=3,
        help='''Relaxation scheme. 3 is default''')


email='nicholaswagner2018@u.northwestern.edu'

args = parser.parse_args()
isif = int(args.isif)
params = args.arguments.split("_")
if params[0] == "b1029":
    ppn = 24
else:
    ppn=20

# Store ENCUT value for changing later
with open('INCAR.static', 'r') as inF:
    for line in inF:
        if 'ENCUT =' in line:
            encut = line.split()[2]

#MOAB commands are always there
with open('run_vasp.q', 'w') as f:
    f.write('#!/bin/bash\n')
    f.write('#MSUB -l nodes=%s:ppn=%s\n' % (params[1], ppn))
    f.write('#MSUB -A %s\n' % params[0])
    f.write('#MSUB -q buyin\n')
    f.write('#MSUB -N %s\n' % params[2])
    f.write('#MSUB -m ae\n')
    f.write('#MSUB -M %s\n' % email)
    f.write('#MSUB -l walltime=%s:00:00\n' % params[3])
    f.write('''cd $PBS_O_WORKDIR
module load mpi/openmpi-1.6.3-intel2013.2
ulimit -s unlimited
nprocs=`wc -l $PBS_NODEFILE | awk '{ print $1 }'`
echo $PBS_NODEFILE\n\n''')

# Boilerplate for DOS and BS calculations
static_text = '''cp INCAR.static INCAR
cp KPOINTS.static KPOINTS
sed -i "s/LREAL = .*/LREAL = .FALSE./" INCAR*
mpirun -n $nprocs /projects/b1027/VASPwannier90.5.3.5/vasp.5.3/vasp > out.static
mv OUTCAR OUTCAR.static
mv vasprun.xml %s_dos.xml
        
sed "s/ICHARG = .*/ICHARG = 11/" INCAR.static > INCAR.bands
sed -i "s/ISMEAR = .*/ISMEAR = 0/" INCAR.bands
sed -i "s/.*SIGMA = .*/SIGMA = 0.1/" INCAR.bands
sed -i "s/.*LCHARG = .*/LCHARG = .FALSE./" INCAR.bands

aflow --kpath < POSCAR > KPOINTS.bands
sed -i "1,/KPOINTS/d" KPOINTS.bands
sed -i '$d' KPOINTS.bands
cp KPOINTS.bands KPOINTS
cp INCAR.bands INCAR
mpirun -n $nprocs /projects/b1027/VASPwannier90.5.3.5/vasp.5.3/vasp > out.bands
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
sed -i "s/.*ISIF = .*/ISIF = %s/" INCAR.is%s.ib1
sed -i "s/.*IBRION = .*/IBRION = 1/" INCAR.is%s.ib1
sed -i "s/.*POTIM = .*/POTIM = 0.2/" INCAR.is%s.ib1
sed -i "s/.*EDIFFG = .*/EDIFFG = -0.0005/" INCAR.is%s.ib1
sed -i "s/ENCUT = .*/ENCUT = %s/" INCAR.is%s.ib1
sed "s/IBRION = 1/IBRION = 2/" INCAR.is%s.ib1 > INCAR.is%s.ib2
sed -i "s/LREAL = .*/LREAL = .FALSE./" INCAR*

for it in 1 2 3
do
cp INCAR.is%s.ib2 INCAR
mpirun -n $nprocs /projects/b1027/VASPwannier90.5.3.5/vasp.5.3/vasp > out.is%s.ib2.$it
mv OUTCAR OUTCAR.is%s.ib2.$it
cp CONTCAR CONTCAR.is%s.ib2.$it
cp  CONTCAR.is%s.ib2.$it POSCAR

cp INCAR.is%s.ib1 INCAR
mpirun -n $nprocs /projects/b1027/VASPwannier90.5.3.5/vasp.5.3/vasp > out.is%s.ib1.$it
mv OUTCAR OUTCAR.is%s.ib1.$it
cp CONTCAR CONTCAR.is%s.ib1.$it
cp CONTCAR POSCAR

done\n\n'''

#Commands for calculating DOS and band structure
if args.jtype == 'static':
    with open('run_vasp.q', 'a') as f:
        f.write(static_text % (params[2], params[2]))


# Commands for full ionic relaxation (volume, positions, and cell parameters)
if args.jtype == 'relax':
    rel_encut = str(int(encut)*0.9)
    with open('run_vasp.q', 'a') as f:
        # Format relax_text with isif number and encut number
        a = [isif]*9
        a.append(rel_encut)
        a.extend([isif]*12)
        f.write(relax_text % tuple(a))
        f.write(static_text % (params[2], params[2]))










