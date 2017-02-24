#!/home/nicholas/CODE/git/aiida/bin/python
import pymatgen as mg
from pymatgen.analysis.elasticity.strain import DeformedStructureSet
import os

structure = mg.Structure.from_file("POSCAR")
def_set = DeformedStructureSet(structure)

count = 1
for struct in def_set:
    struct.to(fmt="POSCAR", filename="./POSCAR%s" % count)
    count+=1
