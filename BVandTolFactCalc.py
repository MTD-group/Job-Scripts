
# coding: utf-8

# In[46]:


# A script to calculate tolerance factors of ABX3 perovskites using bond valences from 2016
# Data from the International Union of Crystallography
# Author: Nick Wagner
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import pymatgen as mg
from pymatgen.analysis.bond_valence import calculate_bv_sum, calculate_bv_sum_unordered, BVAnalyzer
import sys
import re


# In[2]:


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


# In[6]:


# Test using BaMnO3
# Should return 1.09630165911 for perovskite and 1.07615743313 for rp=2
print(calc_tol_factor(['Ba', 'Mn','O'], ['2+', '4+', '2-']))
print(calc_tol_factor(['Ba', 'Mn','O'], ['2+', '4+', '2-'], rp=2))


# In[5]:


def compute_bv(params, pair):
    bv = np.exp((params['Ro']- pair[1])/params['B'])
    return bv


# In[12]:


def isanion(atom, anions=['O', 'S', 'F', 'Cl']):
    #print "in isanion fun... atom is {} and anions are {}".format(atom, anions)
    check = atom in anions
    return check

def iscation(atom, cations=[]):
    check = atom not in ['O', 'S', 'F', 'Cl'] 
    return check    


def gii_compute(struct, formal_val={}):
    '''
    Will compite the global instability index of a structure provided you give a dictionary of valences
    for each species.
    
    This code will not work if you have true non-integer valences.
    
    TODO:
    1. Add support for sites with multiple occupancy but species having only one integer valence
    2. Add ability to specify valence for each site to deal with charge ordering.
    '''
    el = struct.species_and_occu[0].elements[0].symbol
    cutoff = 6
    pymat_neighbors = struct.get_all_neighbors(cutoff, include_index=True)
    if not formal_val:
        print("Please specify formal valences of all species. Returning None")
        return
    
    # for loop to calculate the BV sum on each atom
    bv_diffs = []
    for atom_indx, neighbors in enumerate(pymat_neighbors):
        bv = 0
        for pair in neighbors:
            atom = struct.species[atom_indx].symbol
            neighbor = struct.species[pair[2]].symbol
            
            try:
                if iscation(atom) and isanion(neighbor):
                    params = get_bv_params(cation=atom, anion=neighbor, 
                                           cat_val=float(formal_val[atom]), an_val=float(formal_val[neighbor]))
                    if compute_bv(params, pair) > 0.03 and atom_indx == 0:
                        print(compute_bv(params, pair))
                    bv += compute_bv(params, pair)
                elif iscation(neighbor) and isanion(atom):
                    params = get_bv_params(cation=neighbor, anion=atom, cat_val=float(formal_val[neighbor]), 
                                           an_val=float(formal_val[atom]))
                    bv += compute_bv(params, pair)
            except:
                print("Trouble with atom: {} and neighbor: {} in {}".format(atom, neighbor, 
                                                                            struct.formula))
                print('Looking for +/- 1 similar cation valence states in BV table')
                try:
                    if iscation(atom) and isanion(neighbor):
                        params = get_bv_params(cation=atom, anion=neighbor, 
                                   cat_val=float(formal_val[atom])+1, an_val=float(formal_val[neighbor]))
                        bv += compute_bv(params, pair)
                    elif iscation(neighbor) and isanion(atom):
                        params = get_bv_params(cation=neighbor, anion=atom, 
                                   cat_val=float(formal_val[neighbor])+1, an_val=float(formal_val[atom]))
                        bv += compute_bv(params, pair)
                    formal_val[atom] = float(formal_val[atom])+1
                except:
                    try:
                        if iscation(atom) and isanion(neighbor):
                            params = get_bv_params(cation=atom, anion=neighbor, 
                                   cat_val=float(formal_val[atom])+1, an_val=float(formal_val[neighbor]))
                            bv += compute_bv(params, pair)
                        elif iscation(neighbor) and isanion(atom):
                            params = get_bv_params(cation=neighbor, anion=atom, 
                                   cat_val=float(formal_val[neighbor])+1, an_val=float(formal_val[atom]))
                            bv += compute_bv(params, pair)
                        formal_val[atom] = float(formal_val[atom])-1
                    except:
                        print("No similar valence states found. Returning None")
                        return None
#         print('Atom: {}, BV: {}'.format(struct.species[atom_indx].symbol, bv))

        bv_diffs.append(np.power(abs(float(formal_val[struct.species[atom_indx].symbol])) - bv, 2))
#         print('BV_diffs: {}'.format(bv_diffs))
        
    GII_val = np.sqrt(np.sum(bv_diffs)/struct.composition.num_atoms)
    return GII_val


# In[83]:


def isanion(atom, anions=['O', 'S', 'F', 'Cl']):
    #print "in isanion fun... atom is {} and anions are {}".format(atom, anions)
    check = atom in anions
    return check

def iscation(atom, cations=[]):
    check = atom not in ['O', 'S', 'F', 'Cl'] 
    return check    


def gii_compute(struct, formal_val={}):
    '''
    Will compite the global instability index of a structure provided you give a dictionary of valences
    for each species.
    
    This code will not work if you have true non-integer valences.
    
    Inputs: 
    
    struct     - a Pymatgen structure object
    formal_val - a dictionary mapping each element to its formal valence number
    
    TODO:
    1. Add support for sites with multiple occupancy but species having only one integer valence
    2. Add ability to specify valence for each site to deal with charge ordering.
    '''
    el = struct.species_and_occu[0].elements[0].symbol
    cutoff = 6
    pymat_neighbors = struct.get_all_neighbors(cutoff, include_index=True)
    if not formal_val:
        print("Please specify formal valences of all species. Returning None")
        return
    
    # for loop to calculate the BV sum on each atom
    bv_diffs = []
    for atom_indx, neighbors in enumerate(pymat_neighbors):
        bv = {} # Dictionary for bond valence sum of each species
        atom = struct.species_and_occu[atom_indx].as_dict() # Dictionary of species and 
                                                            # their occupancies at the site
        for key in atom.keys():
            new_name = re.sub(r'[^a-zA-Z ]+', '', key) # Strip out numbers and symbols from the element names
            atom[new_name] = atom.pop(key)
            bv[new_name] = 0
        
        
        for pair in neighbors:
            neighbor = struct.species_and_occu[pair[2]].as_dict()
            for key in neighbor.keys(): 
                new_name = re.sub(r'[^a-zA-Z ]+', '', key)
                neighbor[new_name] = neighbor.pop(key)
                
            for species1 in list(atom): # Iterate over dictionary keys aka elemental species
                for species2 in list(neighbor):
                    try:
                        if iscation(species1) and isanion(species2):
                            params = get_bv_params(cation=species1, anion=species2, 
                                                   cat_val=float(formal_val[species1]), 
                                                   an_val=float(formal_val[species2]))
                            bv[species1] += compute_bv(params, pair)
                        elif iscation(species2) and isanion(species1):
                            params = get_bv_params(cation=species2, anion=species1, 
                                                   cat_val=float(formal_val[species2]), 
                                                   an_val=float(formal_val[species1]))
                            bv[species1] += neighbor[species2]*compute_bv(params, pair)
                    except:
                        print("Trouble with atom: {} and neighbor: {} in {}".format(species1, species2, 
                                                                                    struct.formula))
                        print('Looking for +/- 1 similar cation valence states in BV table')
                        try:
                            if iscation(species1) and isanion(species2):
                                params = get_bv_params(cation=species1, anion=species2, 
                                           cat_val=float(formal_val[species1])+1, 
                                                       an_val=float(formal_val[species2]))
                                bv[species1] += compute_bv(params, pair)
                            elif iscation(species2) and isanion(species1):
                                params = get_bv_params(cation=species2, anion=species1, 
                                           cat_val=float(formal_val[species2])+1, an_val=float(formal_val[species1]))
                                bv[species1] += compute_bv(params, pair)
                            formal_val[species1] = float(formal_val[species1])+1
                        except:
                            try:
                                if iscation(species1) and isanion(species2):
                                    params = get_bv_params(cation=species1, anion=species2, 
                                           cat_val=float(formal_val[species1])+1, an_val=float(formal_val[species2]))
                                    bv[species1] += compute_bv(params, pair)
                                elif iscation(species2) and isanion(species1):
                                    params = get_bv_params(cation=species2, anion=species1, 
                                           cat_val=float(formal_val[species2])+1, 
                                                           an_val=float(formal_val[species1]))
                                    bv[species1] += compute_bv(params, pair)
                                formal_val[species1] = float(formal_val[species1])-1
                            except:
                                print("No similar valence states found. Returning None")
                                return None
                            
                            
                            
                            
#         print('Atom: {}, BV: {}'.format(struct.species[atom_indx].symbol, bv))
        
        weighted_val_diff = 0.0
        for key in list(atom): # Compute the weighted average formal valence in case of multiple occupancy
            formal_valence = abs(float(formal_val[key]))
            weighted_val_diff += atom[key] * (formal_valence - bv[key])
            print("Bond valence of {}: {}".format(key,bv[key]))
        bv_diffs.append(np.power(weighted_val_diff, 2))
        
#         print('BV_diffs: {}'.format(bv_diffs))
        
    GII_val = np.sqrt(np.sum(bv_diffs)/struct.composition.num_atoms)
    return GII_val


# In[67]:


def bv_sum(atom, ):
    if iscation(atom) and isanion(neighbor):
        params = get_bv_params(cation=atom, anion=neighbor, cat_val=float(formal_val[atom]), 
                               an_val=float(formal_val[neighbor]))
        bv += compute_bv(params, pair)
    elif iscation(neighbor) and isanion(atom):
        params = get_bv_params(cation=neighbor, anion=atom, cat_val=float(formal_val[neighbor]), 
                               an_val=float(formal_val[atom]))
        bv += compute_bv(params, pair)
    return bv


# In[69]:


def get_coordination_number(struct, formal_val={}):
    print('This is {}'.format(struct.composition.reduced_formula))
    cutoff = 6
    elem_coord_dict = {}
    for indx, site in enumerate(struct.sites):
        neighbors = struct.get_neighbors(site, cutoff, include_index=True)
        coordnum = 0
        for buddy in neighbors:
            atom = struct.species[indx].symbol
            neighbor = struct.species[buddy[2]].symbol
            if iscation(atom) and isanion(neighbor):
                params = get_bv_params(cation=atom, anion=neighbor, 
                                       cat_val=float(formal_val[atom]), an_val=float(formal_val[neighbor]))
                if compute_bv(params, buddy) > 0.09:
                    coordnum += 1
#                     print(compute_bv(params, buddy))
            elif iscation(neighbor) and isanion(atom):
                params = get_bv_params(cation=neighbor, anion=atom, 
                                       cat_val=float(formal_val[neighbor]), an_val=float(formal_val[atom]))
                if compute_bv(params, buddy) > 0.09:
                    coordnum += 1
        if (site.species_and_occu.alphabetical_formula, coordnum) not in elem_coord_dict.items():            
            print("Index: {}, Species: {}, Coordination: {}".format(indx, 
                                                    site.species_and_occu.alphabetical_formula, coordnum))
            elem_coord_dict[site.species_and_occu.alphabetical_formula] = coordnum


# In[88]:


# Experiment with finding coordination numbers of specific compounds
# BaTiO3 in space group 123 should return a value around 0.367
df = pd.read_excel("../data/IMT_Classification_Dataset.xlsx",sheetname="Combined_MIT+nonMIT")
df = df.loc[df['Compound'].str.len() < 10]
df = df.loc[df['Compound'] == 'LaSrCrO4']
valences = df['formal_val'].str.split("_")
df['unrounded_vals'] = valences.apply(pd.to_numeric)
df['elems'] = df['Elements'].str.split("_")
def get_formal_vals(row):
    return dict(zip(row['elems'], row['rounded_vals'])) 
    
# df['formal_vals'] = df.apply(get_formal_vals, axis=1)

gii_values = []
for i in df.index:
    struct = mg.Structure.from_file('..' + df.loc[i, "struct_file_path"])
    rounded_vals = np.around(df.loc[i, 'unrounded_vals'])
    formal_vals = dict(zip(df.loc[i, 'elems'], rounded_vals))
    gii = gii_compute(struct, formal_vals)
    print('GII: {}'.format(gii))
#     formal_val = df.loc[i, 'formal_vals']
#     try:
#         get_coordination_number(struct, formal_vals)
#     except:
#         print('{} failed. Please manually inspect'.format(df.loc[i, 'Compound']))


# In[18]:


# Compute GII of one cif file
struct = mg.Structure.from_file("../../../../Downloads/EntryWithCollCode173697.cif")
formal_vals = {'Sr': 2, 'Co': 3, 'O': -2}
gii_compute(struct, formal_vals)


# In[19]:


df = pd.read_excel("../data/IMT_Classification_Dataset.xlsx",sheetname="Combined_MIT+nonMIT")
df.gii.describe()


# In[10]:


struct.composition.reduced_formula


# In[ ]:


# Experiment with specific compounds
# BaTiO3 in space group 123 should return a value around 0.367
df = pd.read_excel("../data/Dataset.xlsx",sheetname="Combined_MIT+nonMIT")
df = df.loc[df['Compound'].str.len() == 5]
valences = df['formal_val'].str.split("_")
df['rounded_vals'] = valences.apply(pd.to_numeric)
df['elems'] = df['Elements'].str.split("_")
def get_formal_vals(row):
    return dict(zip(row['elems'], row['rounded_vals'])) 
    
df['formal_vals'] = df.apply(get_formal_vals, axis=1)

gii_values = []
for i in df.index:
    struct = mg.Structure.from_file('..' + df.loc[i, "struct_file_path"])
    gii = gii_compute(struct, formal_val)
    gii_values.append(gii)
    formal_val = df.loc[i, 'formal_vals']
    get_coordination_number(struct, formal_vals)
    
foo = pd.DataFrame(gii_values, columns=['GII'])
foo['Compound'] = df['Compound'].values
foo


# In[16]:


# Calculate GII for all compounds in your dataset
# Requires one column with elements (e.g. Ba_Ti_O),
# one column with the structure path (e.g. ./Structures/BaTiO3.cif),
# and one column with the formal valences (e.g. 2_4_-2).
# Does not work with disordered compounds due to a Pymatgen limitation
# Sites with mixed valence will have their valence state rounded to the nearest integer
# Output saved to GII_temp.csv in data folder
try:
    df = pd.read_excel("../data/Dataset.xlsx",sheetname="Combined_MIT+nonMIT", 
                       usecols=['Compound', 'formal_val', 'Elements', 'struct_file_path'])
#     df = df.loc[df['Insulator'] < 2]
except:
    print("Define df to be your dataset path in the code!")

gii_values = []
for i in df.index:
    struct = mg.Structure.from_file('..' + df.loc[i, "struct_file_path"])
    valences = df.iloc[[i]]['formal_val'].tolist()[0].split("_")
    rounded_vals = np.around(np.array(valences).astype(float))
    formal_val = dict(zip(df.iloc[[i]]['Elements'].tolist()[0].split("_"), 
                      rounded_vals))
    name=df.loc[i,'Compound']
    try:
        gii = gii_compute(struct, formal_val)
        gii_values.append(gii)
    except AttributeError:
        print("Attribute Error. You might have mixed occupations here in: {}".format(name))
        gii_values.append(None)
df['gii'] = gii_values
df.to_csv("../data/GII_temp.csv")    


# In[ ]:


# Generate Prasanna's chemistries from doi:10.1038/ncomms14282 and calculate their tolerance factors
# for Danilo
chem_ids = pd.read_excel("../../RP_tolerance_factors/Site_IDs.xlsx", sheetname='Sheet1')
asites = chem_ids.ix[:, 0:3]
bsites = chem_ids.ix[0:25, 3:]
# print(asites)
def generate_chem(asites, bsites):
    import itertools
    acombos = list(itertools.combinations(asites.AsiteID, 2)) # Find unique combinations of A and A' sites
    df = pd.DataFrame(columns=['A1', 'A2', 'B', 'X', 'A1_valence', 'A2_valence', 'B_valence', 'X_valence'])
    for apair in acombos:        
        a1_info = asites[asites['AsiteID'] == apair[0]][['ElemA', 'Avalence']]
        a2_info = asites[asites['AsiteID'] == apair[1]][['ElemA', 'Avalence']]
        
        a1_info.rename(index=str, columns={'ElemA': 'A1', 'Avalence': 'A1_valence'}, inplace=True)
        a2_info.rename(index=str, columns={'ElemA': 'A2', 'Avalence': 'A2_valence'}, inplace=True)
        a_sites = a1_info
        a_sites['A2'] = a2_info['A2'].values
        a_sites['A2_valence'] = a2_info['A2_valence'].values
        for index, bsite in bsites.iterrows():
            row = a_sites
            row['B'] = bsite['ElemB']
            row['B_valence'] = bsite['Bvalence']
            row['X'] = 'O'
            row['X_valence'] = -2
            df = df.append(row) 
    return df        
#                       row.ElemB, row.Bvalence

            
crazy_chems = generate_chem(asites, bsites)
sane_chems = crazy_chems[crazy_chems['A1_valence'] + crazy_chems['A2_valence'] + crazy_chems['B_valence'] == 8]
sane_chems = sane_chems[(sane_chems['A1'] != 'Pm') & (sane_chems['A2'] != 'Pm')]
sane_chems.to_csv("./Compounds_forDanilo.csv")
sane_chems.tail()    


# In[ ]:


# Generate Prasanna's chemistries from doi:10.1038/ncomms14282 for the Mobilities Project for Ken's perusal
chem_ids = pd.read_excel("../../RP_tolerance_factors/Site_IDs.xlsx", sheetname='Sheet2') # Sheet 2 is reduced
if chem_ids.shape[0] == 30:
    asites = chem_ids[['ElemA', 'Avalence']].iloc[:]
    bsites = chem_ids[['ElemB', 'Bvalence']].iloc[0:26]
else:
    asites = chem_ids[['ElemA', 'Avalence']].iloc[0:14]
    bsites = chem_ids[['ElemB', 'Bvalence']].iloc[:]


def generate_perov_chem(asites, bsites):
    import itertools
    df = pd.DataFrame(columns=['A', 'B', 'X', 'A_valence', 'B_valence', 'X_valence'])        
    asites.rename(index=str, columns={'ElemA': 'A', 'Avalence': 'A_valence'}, inplace=True)
    for index, asite in asites.iterrows():

        for index, bsite in bsites.iterrows():
            row = asite
            row['B'] = bsite['ElemB']
            row['B_valence'] = bsite['Bvalence']
            row['X'] = 'O'
            row['X_valence'] = -2
            df = df.append(row) 
    return df        

            
crazy_chems = generate_perov_chem(asites, bsites)
sane_chems = crazy_chems[crazy_chems['A_valence'] + crazy_chems['B_valence'] == 6]
sane_chems = sane_chems[(sane_chems['A'] != 'Pm')]
if chem_ids.shape[0] == 30:
    sane_chems.to_csv("../data/PotentialCompounds_extended.csv")
else:
    sane_chems.to_csv("../data/PotentialCompounds.csv")
sane_chems.tail()    


# In[ ]:





# In[ ]:


sane_chems.head()


# In[ ]:


# Calculate tolerance factors for Danilo
sane_chems = pd.read_csv("./Compounds_forDanilo.csv", index_col=0)
df = sane_chems

tol_factors = []
for i in range(len(sane_chems.index)):
    ion_list = [df.loc[i, 'A1'], df.loc[i, 'A2'], df.loc[i, 'B'], df.loc[i, 'X']]
    valence_list = [int(df.loc[i, 'A1_valence']), int(df.loc[i, 'A2_valence']), 
                    int(df.loc[i, 'B_valence']), int(df.loc[i, 'X_valence'])]
    try:        
        tol_fact = calc_tol_factor(ion_list, valence_list, rp=1)
    except(AttributeError):
        print("Your compound: {} had an issue and will not be calculated".format(ion_list))
    tol_factors.append(tol_fact)
    
df_tol = pd.DataFrame(tol_factors, columns=['Tol_factor'])

df['Tol_factor'] = df_tol
df.head()


# In[ ]:


df.to_csv("./Compounds_forDanilo.csv")


# In[ ]:


df.describe()


# In[ ]:


# df = pd.DataFrame(columns=['A1', 'A2', 'B', 'X', 'A1_valence', 'A2_valence', 'B_valence', 'X_valence'])
# arr = np.array([10, 20]).reshape(1,2)
# foo = pd.DataFrame(arr, columns=['A1','A1_valence'])
# df = df.append(foo)
df.tail()


# In[ ]:


# Test using BaMnO3
print(calc_tol_factor(['Ba', 'Mn','O'], ['2+', '4+', '2-']))
print(calc_tol_factor(['Ba', 'Mn','O'], ['2+', '4+', '2-'], rp=2))


# In[ ]:


# Compare vanadate tolerance factors calculated with my code to that of Nicole Benedek's
names = ['Yb','Dy','Ho','Y','Tb','Gd', 'Eu','Sm','Nd','Pr','Ce','La']
nicole = [0.883, 0.901, 0.897, 0.827, 0.906, 0.912, 0.916, 0.92, 0.93, 0.936, 0.942, 0.95]
nick = []
for name in names:
    nick.append(float('{:0.3f}'.format(calc_tol_factor([name, 'V','O'], ['3+', '3+', '2-']))))
d = {'nicole': nicole, 'nick': nick}
vanadates = pd.DataFrame(data=d ,index=names)
vanadates


# In[ ]:


ax = vanadates.plot.bar(figsize=(16,14),fontsize=32)
ax.set_ylabel('Tolerance Factor', fontsize=32)
ax.set_ylim(0.8, 1)
ax.legend(fontsize=32)
ax.set_title('Vanadate Tolerance Factors', fontsize=36)
plt.show()


# In[ ]:


nickels = ['Lu', 'Y', 'Dy', 'Gd', 'Eu', 'Sm', 'Nd', 'Pr', 'La']
nicole = [0.904, 0.851, 0.928, 0.938, 0.942, 0.947, 0.957, 0.964, 0.977]
nick= []
for nickel in nickels:
    nick.append(float('{:0.3f}'.format(calc_tol_factor([nickel, 'Ni','O'], ['3+', '3+', '2-']))))

d = {'nicole': nicole, 'nick': nick}
nickelates = pd.DataFrame(data=d, index=nickels)
ax = nickelates.plot.bar(figsize=(16,14),fontsize=32)
ax.set_ylabel('Tolerance Factor', fontsize=32)
ax.set_ylim(0.8, 1)
ax.legend(fontsize=32)
ax.set_title('Nickelate Tolerance Factors', fontsize=36)
plt.show()

