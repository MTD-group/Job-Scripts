
# coding: utf-8

# In[2]:

# A script to calculate tolerance factors of ABX3 perovskites using bond valences from 2016
# Data from the International Union of Crystallography
# Author: Nick Wagner
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import pymatgen as mg
from pymatgen.analysis.bond_valence import calculate_bv_sum, calculate_bv_sum_unordered, BVAnalyzer


# In[ ]:

bv = pd.read_csv("../data/Bond_valences2016.csv")
bv.head()


# In[3]:

# Use element names and valences to lookup bond valence
def get_bv_params(cation, anion, cat_val, an_val):
    bond_val_list = bv[(bv['Atom1'] == cation) & (bv['Atom1_valence'] == cat_val)                  & (bv['Atom2'] == anion) & (bv['Atom2_valence'] == an_val)]
    return bond_val_list.iloc[0] # If multiple values exist, take first one


# In[4]:

# A function to calculate a generalized Goldschmidt tolerance factor for perovskites and RP phases
def calc_tol_factor(ion_list, valence_list, rp=0):
    if len(ion_list) > 4 or len(ion_list) < 3:
        print("Error: there should be three or four elements")
        return None
    if len(ion_list) < 4:
        for i in range(len(valence_list)): # If charge is 2-, make -2 to match tables
            if valence_list[i][-1] == '-':
                valence_list[i] = valence_list[i][-1] + valence_list[i][:-1]
        for i in range(len(valence_list)): # Similarly, change 2+ to 2
            valence_list[i] = int(valence_list[i].strip("+"))
        
    if len(ion_list) == 4:
#         print("RED ALERT: We are taking averages of bond valence parameters")
        AO_value1 = get_bv_params(ion_list[0], ion_list[-1], valence_list[0], valence_list[-1])
        AO_value2 = get_bv_params(ion_list[1], ion_list[-1], valence_list[1], valence_list[-1])
        AO_values = np.concatenate([AO_value1.values.reshape(1, len(AO_value1)), 
                                    AO_value2.values.reshape(1, len(AO_value2))])
        AO_B = np.average(AO_values[:, 4])
        AO_Ro = np.average(AO_values[:, 5])
        AO_valence = np.average(AO_values[:, 1]) # RED ALERT: We are taking averages of bond valence parameters
    else:
        AO_row = get_bv_params(ion_list[0], ion_list[-1], valence_list[0], valence_list[-1])
    
    BO_row = get_bv_params(ion_list[-2], ion_list[-1], valence_list[-2], valence_list[-1])
    
    
    if len(ion_list) != 4:
        if rp == 0:
            AO_bv = AO_row['Ro']-AO_row['B'] * np.log(AO_row['Atom1_valence']/12)
            BO_bv = BO_row['Ro']-BO_row['B'] * np.log(BO_row['Atom1_valence']/6)               
        else: # Currently for Ruddlesden-Popper phases a naive weighted sum is used between A-site coordination of 
              # 9 in the rocksalt layer and 12 in perovskite
            AO_bv = AO_row['Ro']-AO_row['B'] * np.log(AO_row['Atom1_valence']/((9+12*(rp-1))/rp))
            BO_bv = BO_row['Ro']-BO_row['B'] * np.log(BO_row['Atom1_valence']/6)
    else:
        if rp == 0:
            AO_bv = AO_Ro-AO_B * np.log(AO_valence/12)
            BO_bv = BO_row['Ro']-BO_row['B'] * np.log(BO_row['Atom1_valence']/6)               
        else: # Currently for Ruddlesden-Popper phases a naive weighted sum is used between A-site coordination of 
              # 9 in the rocksalt layer and 12 in perovskite
            AO_bv = AO_Ro-AO_B * np.log(AO_valence/((9+12*(rp-1))/rp))
            BO_bv = BO_row['Ro']-BO_row['B'] * np.log(BO_row['Atom1_valence']/6)
    
    tol_fact = AO_bv / (2**0.5 * BO_bv)
    
    return tol_fact


# In[ ]:

# Test using BaMnO3
# Should return 1.09630165911 for perovskite and 1.07615743313 for rp=2
print(calc_tol_factor(['Ba', 'Mn','O'], ['2+', '4+', '2-']))
print(calc_tol_factor(['Ba', 'Mn','O'], ['2+', '4+', '2-'], rp=2))


# In[5]:

def isanion(atom, anions=['O', 'S', 'F', 'Cl']):
    #print "in isanion fun... atom is {} and anions are {}".format(atom, anions)
    check = atom in anions
    return check

def iscation(atom, cations):
    check = atom not in ['O', 'S', 'F', 'Cl'] 
    return check    


def MObonds_greedy(structure,Msite, cutoff=3.0):
# This function takes a pymatgen structure and perovskite Bsite and returns a list of the bond lengths associated with the Msite/anion bond lengths for the first site
    bond_lengths = []
    # determine Bsite and oxygen indexes
    for site in structure.sites:
        if Msite in str(site):
            neighbors = structure.get_neighbors(site, r = cutoff, include_index=True)
            for neighbor in neighbors:
                elems_on_neighsite = structure.species_and_occu[neighbor[2]].elements
                symbols = [elem.symbol for elem in elems_on_neighsite]
                if Msite in symbols:
                    continue
                else:
                    bond_lengths.append(neighbor[1])
            if not bond_lengths:
                neighbors = structure.get_neighbors(site, r = cutoff+0.6, include_index=True)
                for neighbor in neighbors:
                    elems_on_neighsite = structure.species_and_occu[neighbor[2]].elements
                    symbols = [elem.symbol for elem in elems_on_neighsite]
                    if Msite in symbols:
                        continue
                    else:
                        bond_lengths.append(neighbor[1])
                return bond_lengths
            else:
                return bond_lengths
    return bond_lengths

# Computes GII using a automatic cutoff determining scheme. The cutoff is intended to be the longest nearest-neighbor bond length
def gii_compute(struct, name):
    el = struct.species_and_occu[0].elements[0].symbol
    max_MObond = np.max(MObonds_greedy(struct, el))
    cutoff = max_MObond
    
    # for loop to calculate the BV sum on each atom
    BVpara = pd.read_csv("../data/Bond_valences2016.csv")
    bv = BVAnalyzer(max_radius=cutoff+0.1)
    bv_diffs = []
    for atom_indx, site in enumerate(struct.sites):
        neighbors = struct.get_neighbors(site, cutoff)
        try:
            bv_sum = calculate_bv_sum_unordered(site, neighbors)
        except:
            bv_sum = calculate_bv_sum(site, neighbors)
        try:
            formal_val = bv.get_valences(struct)[atom_indx]
        except:
            print('Difficulty obtaining site valences. Returning None')
            return None
        bv_diffs.append(np.power(np.subtract(bv_sum, formal_val),2))
    GII_val = np.sqrt(np.sum(bv_diffs)/struct.composition.num_atoms)
    return GII_val


# In[6]:

# Calculate GII for all compounds in your dataset
# Requires one column with elements (e.g. Ba_Ti_O),
# one column with the structure path (e.g. ./Structures/BaTiO3.cif),
# and one column with the formal valences (e.g. 2_4_-2).
# Does not work with disordered compounds due to a Pymatgen limitation
# Output saved to GII_temp.csv in current folder.
# Ideally GII values should be close to 0 and less than 0.2, but good luck
try:
    df = pd.read_excel("../data/Dataset.xlsx",sheetname="Combined_MIT+nonMIT")
except:
    print("Define df to be your dataset path in the code!")
gii_values = []
for i in df.index:
    struct = mg.Structure.from_file('..' + df.loc[i, "struct_file_path"])
    name=df.loc[i,'Compound']
    gii = gii_compute(struct, name)
    gii_values.append(gii)
    
foo = pd.DataFrame(gii_values)
foo.to_csv("../data/GII_temp.csv")    


