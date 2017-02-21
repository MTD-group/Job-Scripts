# Call with your python executable to generate combined DOS/bandstructure figure 
# from your two vasprun files
# Relies upon pymatgen version 4.5.6 or greater, matplotlib, and glob
# Also outputs the size and type of bandgap to 'bandgap.txt' and prints this to std out

from pymatgen.io.vasp.outputs import Vasprun
from pymatgen.electronic_structure.plotter import BSDOSPlotter
import glob
import matplotlib as plt

bs_name = sorted(glob.glob('*.xml'))[0] # Takes the xml files in alphabetical order  
dos_name = sorted(glob.glob('*.xml'))[1] # e.g. bs.xml, dos.xml


bs_vrun = Vasprun(bs_name, parse_projected_eigen=True)
dos_vrun = Vasprun(dos_name)
temp = bs_vrun.get_band_structure(kpoints_filename = 'KPOINTS.bands')
vbm = temp.get_vbm()['energy'] # Get valence band max for better looking plots of semiconductors

bs = bs_vrun.get_band_structure(kpoints_filename = 'KPOINTS.bands', efermi=vbm)
dos = dos_vrun.complete_dos

with open('bandgap.txt', 'w') as f:
    f.write('Bandgap: ' + str(bs.get_band_gap()))
print('Bandgap: ' + str(bs.get_band_gap()))
bsdosplot = BSDOSPlotter()
#If compound is metallic, leave fermi level at what VASP says, otherwise set zero at valence band max
if temp.is_metal():
    figuro = bsdosplot.get_plot(temp, dos)
else:
    figuro = bsdosplot.get_plot(bs, dos)
figuro.savefig('bsdos.png')
