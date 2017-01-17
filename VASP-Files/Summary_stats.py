'''
Outputs some important values for a DFT calculation: Magnetization 
of ions, elastic tensor if calculated, and final energy
Author: Nicholas Wagner (MTD Group)
'''
from pymatgen.io.vasp.outputs import Outcar
from pymatgen.io.vasp.outputs import Vasprun
import glob
import pprint

out = Outcar('OUTCAR.static')
mag = out.magnetization

with open('summary.txt', 'w') as f:
    f.write('\n\nMagnetization: \n')
    pp = pprint.PrettyPrinter(indent=4, stream=f)
    pp.pprint(mag)

    try:
        elasts = out.elastic_tensor
        f.write('Elastic tensor: \n')
        pp.pprint(elasts)
    except:
        f.write('No elastic tensor found\n')

    try:
        vrun = Vasprun(sorted(glob.glob('*.xml'))[1], parse_dos=False, parse_eigen=False)
        f.write('Final energy: ' + str(vrun.final_energy)+'\n\n\n')
    except:
        try:
            vrun = Vasprun(sorted(glob.glob('*.xml'))[0], parse_dos=False, parse_eigen=False)
            f.write('Final energy: ' + str(vrun.final_energy)+'\n\n\n')
        except:
            f.write('No final energy found in vasprun.xml\n\n')
with open('summary.txt', 'r') as f:
    print(f.read())
