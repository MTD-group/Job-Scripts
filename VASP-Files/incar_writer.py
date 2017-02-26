import argparse, os.path

# Author: Nick Wagner
# This file creates VASP INCAR.static files for some common cases
# found with transition metal oxides. Only electronic control parameters
# need to be specified here. Ionic relaxation tags are handled by job
# submission script writers.

parser = argparse.ArgumentParser()
parser.add_argument("-m", action='store', dest='mtype', default='non_mag',
        help='''Magnetism type. "non_mag" and "ferro" are current options. Default
        is "non_mag"''')
parser.add_argument("-name", action='store', dest='system_name', default='noname',
        help='''System name. Default is "noname"''')
args = parser.parse_args()


# Text for a simple non-magnetic calculation of DOS
nonmag_text = '''System = {0}
#Start parameters
ISTART =      1  #  job   : 0-new  1-cont  2-samecut
ICHARG =      1  #  charge: 1-file 2-atom 10-const
INIWAV =      1  #  electr: 0-lowe 1-rand

# Convergence tolerance
EDIFF = 1e-7
SYMPREC = 1e-4
#ALGO = FAST
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
ENCUT = 625
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

mag_text = '''#Magnetism
ISPIN = 2
#MAGMOM = 8*0 4.0*8 24*0
'''

if args.mtype == 'non_mag':
    with open('INCAR.static', 'w') as f:
        f.write(nonmag_text.format(args.system_name))

if args.mtype == 'ferro':
    with open('INCAR.static', 'w') as f:
        f.write(nonmag_text.format(args.system_name))
        f.write(mag_text)