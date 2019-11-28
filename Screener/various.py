import pandas as pd
import numpy as np
import phonopy
# from ase import io
import os
from numpy import genfromtxt

def Make_Sample_set_for_MCEScreening_tests():
    df = pd.read_csv('Alldata.csv', index_col=0, sep=',')
    monoclinic = df[df['lattice_system_relax'].str.match('monoclinic')]  # no preset
    triclinic = df[df['lattice_system_relax'].str.match('triclinic')]  # no preset
    rhombohedral = df[df['lattice_system_relax'].str.match('rhombohedral')]  # no preset
    hexagonal = df[df['lattice_system_relax'].str.match('hexagonal')]  # 3  set       Fe2P
    orthorhombic = df[df['lattice_system_relax'].str.match('orthorhombic')]  # 3  set       MnCoP
    cubic = df[df['lattice_system_relax'].str.match('cubic')]  # 3  set       FeRh
    tetragonal = df[df['lattice_system_relax'].str.match('tetragonal')]  # 1  set       LaFe9Si4

    k = tetragonal
    out = k.sample(n=4, random_state=1)
    out.to_csv('100/tetragonal.csv')

def Calculate_mag_field_for_aflow_datalist():
    # I should rewrite it using pandas I guess - now it;s a bit sketchy
    my_data = genfromtxt('Alldata.csv', delimiter=',', fmt=str)
    i = 0
    while i <= 8043:
        my_data[i, 6] = str(((float(my_data[i, 5])*9.2741*10**(-24))/(float(my_data[i, 4])*10**(-30)))*4*3.14*10**(-7))
        i = i + 1
    np.savetxt("new_data.txt", my_data, delimiter=',')


def Calculate_unique_mag_positions_and_add_to_datalist():
    df = pd.read_csv('data.txt', index_col=0, sep=',')
    for item in df.index.tolist():
        ph = phonopy.load(unitcell_filename="./data/"+str(item))
        sym_list = ph.get_symmetry().get_Wyckoff_letters()
        at_type_len = len(ph.unitcell.get_positions())
        at_list = []
        at_type_list = []

        poscar = open("./data/"+str(item), "r")
        list = []
        for line in poscar:
            list.append(line)
        poscar.close()
        for i in range(7, 7 + at_type_len):
            l = str.split(list[i], '  ')
            at_type_list.append(str(l[4]))

        for num, val in enumerate(sym_list):  ### build list of magnetic positions
            if at_type_list[num] in magnetic:
                at_list.append(val)
        num_unique = len(set(at_list))  ### check how many are unique
        df.loc[item,'mag_sites'] = num_unique  ### return resulting number of sites



# for num, val in enumerate(sym_list):
# if at_type_list[num] in Rb:
#     df = df.drop([item], axis=0)


# print(ph.get_symmetry().get_dataset())

# for i, vall in enumerate(at_type_list):
#     sym_point = ph.get_symmetry().get_site_point_group()
#     # sym_list = sym_list.append(ph.get_symmetry().get_pointgroup(sym_point))
#     print(sym_point)
# print(sym_list)
# print(at_type_list)

# read element list - if Rb is there delete entry and file - later compare by hand if number of deleted is correct
# run my site search app feeding files and write resulting number in the pandas csv (maybe just seperate txt to be easier then jus join in excel
# need to remove last two zero colums???
