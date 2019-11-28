import pandas as pd
import numpy as np
import phonopy
from ase import io
import os
from numpy import genfromtxt

from Screener import sym_detector

One need to make it progressivy going and properly separate functions

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

    my_data = genfromtxt('datalist.csv', delimiter=',', fmt=str)
    i = 0
    while i <= 8043:
        my_data[i, 6] = str(((float(my_data[i, 5])*9.2741*10**(-24))/(float(my_data[i, 4])*10**(-30)))*4*3.14*10**(-7))
        i = i + 1
    np.savetxt("new_data.txt", my_data, delimiter=',')

    # this should be a better approach -test
    df = pd.read_csv('datalist.csv', index_col=0, sep=',')
    for item in df.index.tolist():
        unique = sym_detector.mag_sites_calculator(item['mag_field'])
        df.loc[item,'mag_field'] = unique   # write resulting number of sites into datalist


def Calculate_unique_mag_positions_and_add_to_datalist():
    df = pd.read_csv('datalist.csv', index_col=0, sep=',')
    for item in df.index.tolist():
        unique = sym_detector.mag_sites_calculator(item['count'])
        df.loc[item,'mag_sites'] = unique   # write resulting number of sites into datalist


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
