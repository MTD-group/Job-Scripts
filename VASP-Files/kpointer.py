#!/home/naw356/bin/python
# Script to generate KPOINTS file for SCF calculations
# Requires PRECALC and getKPoints files to be located in ~/bin/ 
import argparse, subprocess, shutil, os.path

# Can be called directly to write a VASP automatic scheme KPOINTS file
if __name__ == '__main__':
    parser = argparse.ArgumentParser()
    parser.add_argument("-n", nargs='+', action="store",
        dest="kpoints", default="6_6_6_G",
        help='''Specify three k point subdivisions
        and mesh type separated by underscores: 
        e.g. 6_6_6_G or 4_2_4_M''')
    args = parser.parse_args()
    params = args.kpoints.split('_')
    if params[3].upper() == 'M':
        mesh = 'Monkhorst-Pack mesh'
    else:
        mesh = 'Gamma centered mesh'

    with open('KPOINTS.static', 'w') as f:
        f.write('''Automatic mesh
0          ! Auto k-point generation scheme
%s          ! %s 
%s %s %s      ! Dimensions''' % (params[3], mesh, params[0], params[1],
        params[2]))

# If called from elsewhere, will use the Mueller group's K-Point Grid Server
# Please cite:  Pandu Wisesa, Kyle A. McGill, and Tim Mueller. 
#               Phys. Rev. B 93, 155109 (2016)
else:
    home = os.path.expanduser("~/bin")
    shutil.copyfile(home + "/PRECALC", "./PRECALC")
    subprocess.call(home + "/getKPoints", shell=True)
    shutil.copyfile("./KPOINTS", "./KPOINTS.static")
