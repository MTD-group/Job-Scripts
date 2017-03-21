import argparse, os.path

# Author: Nick Wagner
# This file creates VASP INCAR.static files for some common cases
# found with transition metal oxides. Only electronic control parameters
# need to be specified here. Ionic relaxation tags are handled by job
# submission script writers.

parser = argparse.ArgumentParser()

# User can set system name if they want
with open("POSCAR") as f:
    poscar_name = f.readline().replace(" ","")
parser.add_argument("-name", action='store', dest='system_name', default=poscar_name,

# Let user set magnetization
parser.add_argument("-m", action='store', dest='magmom',
        help='''MAGMOM tag separated by '_': e.g. 2*0_2_-2_6*0; ISPIN is set
        to 1 if not provided.''')

# Hubbard U settings
parser.add_argument("-hu", action='store', dest='hubU', default='False',
        help='''Whether to use Dudarev's method to include onsite energy term.
        Defaults to false.''')

args = parser.parse_args()
magmom = args.magmom.replace("_", " ")
hubU = args.hubU

# Text for a simple non-magnetic calculation of DOS
nonmag_text = '''System = {0}
#Start parameters
ISTART =      1  #  job   : 0-new  1-cont  2-samecut
ICHARG =      1  #  charge: 1-file 2-atom 11-const

# Convergence tolerance
EDIFF = 1e-7
SYMPREC = 1e-3
NELM = 40

# Density of States Info
NEDOS = 3001
EMIN = -12
EMAX = 12
LORBIT = 11

# BZ Integrations
ISMEAR = -5
#SIGMA = 0.1 # Smearing width for Gaussian

# Calculation Details
LWAVE = .FALSE. # Do we want to write the wavefunctions
LCHARG = .TRUE. # Do we want to write the charge density
PREC = HIGH # ENCUT set according to PAW
ENCUT = 600
LREAL = .FALSE.
LASPH = .TRUE.

# Relaxation parameters
NSW = 0 # number of ionic steps
#ISIF = 3 # calculate force and relax ions
#IBRION = 1 # Newtonian method for mini_opt
#POTIM = 0.3 # quasi_newton scalar
#EDIFFG = -0.0005

#Parallelization
LPLANE = .TRUE.
NPAR   = 4   # This should be a factor of # cores approximately equal to (# cores)^1/2
LSCALU = .FALSE.
NSIM   = 4
GGA = PS # Use PBE-sol XC functional\n'''

# Spin polarized parameters
mag_text = '''#Magnetism
ISPIN = 2
#MAGMOM = {0}
'''

# Hubbard U parameters
hubU_text = '''#LSDA + U
LDAU = .TRUE.
LDAUTYPE = 2
LDAUL = {0} 
LDAUU = {1} 
LDAUJ = 0 0 0 
LMAXMIX = {2}
LDAUPRINT = 1

'''
if not magmom:
    with open('INCAR.static', 'w') as f:
        f.write(nonmag_text.format(args.system_name))
else:
    if hubU is False:
        with open('INCAR.static', 'w') as f:
            f.write(nonmag_text.format(args.system_name))
            f.write(mag_text.format(magmom))



