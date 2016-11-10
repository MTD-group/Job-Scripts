#!/home/naw356/bin/python
# Put your Python executable location on the first line
import argparse, os.path

parser = argparse.ArgumentParser()
parser.add_argument("-n", nargs='+', action="store",
        dest="kpoints", default="6_6_6_G",
        help='''Specify three k point subdivisions
        and mesh type separated by underscores: 
        e.g. 6_6_6_G or 4_2_4_M''')
args = parser.parse_args()
params = args.kpoints.split('_')
if params[3] == 'M':
    mesh = 'Monkhorst-Pack mesh'
else:
    mesh = 'Gamma centered mesh'

with open('KPOINTS.static', 'w') as f:
    f.write('''Automatic mesh
0          ! Auto k-point generation scheme
%s          ! %s 
%s %s %s      ! Dimensions''' % (params[3], mesh, params[0], params[1],
    params[2]))


