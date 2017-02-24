import pymatgen as mg
from pymatgen.analysis.elasticity.strain import DeformedStructureSet
import os
from shutil import copyfile
from pymatgen.io.vasp.outputs import Vasprun
from pymatgen.analysis.elasticity.stress import Stress
from pymatgen.analysis.elasticity.elastic import ElasticTensor

structure = mg.Structure.from_file("POSCAR")
def_set = DeformedStructureSet(structure)
strains = def_set.as_strain_dict()

calculations = []
for x in range(1,25):
    calculations.append('poscar%s' % x)

match_dict = {}
for calc in calculations:
    struct = mg.Structure.from_file(calc + "/POSCAR")
    vrun =  Vasprun(calc + '/vasprun.xml', parse_dos=False, parse_eigen=False)
    stress = Stress(vrun.ionic_steps[-1]['stress'])
    for strain in strains:
        if strains[strain].lattice == struct.lattice:
            match_dict[strain] = stress
elastics = ElasticTensor.from_stress_dict(match_dict)
with open("elasts.txt", 'w') as f:
    f.write(str(elastics.voigt) + '\n')
    f.write(str(elastics.k_voigt) + '\n')
