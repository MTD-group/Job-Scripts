# Requires: pymatgen module
# Generates KPOINTS.bands file with high symmetry lines added for
# an arbitrary POSCAR with symmetry precision=1e-3

import pymatgen as mg
from pymatgen.symmetry.bandstructure import HighSymmKpath
from pymatgen.io.vasp.inputs import Kpoints

struct = mg.Structure.from_file("POSCAR")
kpath = HighSymmKpath(struct)
kpoints = Kpoints.automatic_linemode(16, kpath)
kpoints.write_file("KPOINTS.bands")
