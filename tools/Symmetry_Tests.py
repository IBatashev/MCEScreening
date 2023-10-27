import pandas as pd
from pymatgen.io.vasp import Poscar
import os
import tqdm
import warnings


# This module is used to check if we perform the deformations the way we expect and that VASP would recognize symmetries
# before and after deformations the same way. This is necessary, because aflolib has some structure files in
# non-conventional way, which may cause the function we use to create deformations to produce incorrect results.
# Checks are performed in several steps, after which a single .csv file containing all problematic entries is created in
# this file True or False is listed for all deformations to show which match with undeformed symmetries and which do not

def calc_list(calc_dict):
    sub_calculations_list = []
    if 'undeformed' in list(calc_dict):
        sub_calculations_list.append('undeformed')
    if 'Applied_Field' in list(calc_dict):
        sub_calculations_list.append('Applied_Field')
    if 'volumetric' in list(calc_dict):
        for i in calc_dict['volumetric']:
            sub_calculations_list.append('V_' + str(i))
    if 'uniaxial' in list(calc_dict):
        for i in calc_dict['uniaxial']:
            par_list = ['a_', 'b_', 'c_']  # we expect to always have a set of calculations with all symmetry types, even if not, does not matter
            calc_list = [x + i for x in par_list]
            for j in calc_list:
                sub_calculations_list.append(str(j))
    return sub_calculations_list


def deformation_check(datalist, datadir):
    """Function that reads through inputdir and uses pymatgen to determine the deformation size in % compared to undeformed structure
    for all deformations and writes result to a new .csv file.
    Only works if undeformed subfolders are in the inputdir"""
    df = pd.read_csv(datalist, index_col=0, sep=',')
    warnings.filterwarnings("ignore")
    with tqdm.tqdm(total=len(df.index)) as pbar:  # A wrapper that creates nice progress bar
        pbar.set_description("Processing datalist")
        for item in df.index.tolist():
            pbar.update(1)
            try:
                poscar = Poscar.from_file(datadir + str(item) + '/undeformed/' + 'POSCAR')
                structure = poscar.structure
                a_u = structure.lattice.a
                b_u = structure.lattice.b
                c_u = structure.lattice.c
            except:
                print('Cannot find undeformed structure in input folder for ', item)

            for deformation in os.listdir(datadir + str(item)):
                if deformation == 'undeformed':
                    pass
                elif not any(ext in str(deformation) for ext in ['a', 'b', 'c']):
                    pass
                else:
                    try:
                        poscar = Poscar.from_file(datadir + str(item) + '/' + deformation + '/' + 'POSCAR')
                        structure = poscar.structure
                        if 'a' in str(deformation):
                            def_percent = round(100 * (1 - (structure.lattice.a / a_u)), 2)
                        elif 'b' in str(deformation):
                            def_percent = round(100 * (1 - (structure.lattice.b / b_u)), 2)
                        elif 'c' in str(deformation):
                            def_percent = round(100 * (1 - (structure.lattice.c / c_u)), 2)
                        else:
                            pass
                        df.loc[item, str(deformation)] = str(def_percent)
                    except:
                        print("Something is wrong with ", item, ' ', deformation,
                              " and I don't really know why, maybe structure file is not recognized by pymatgen")
                        df.loc[item, str(deformation)] = "pymatgen failure"
    df.to_csv((datalist.replace('.csv', '_def_check.csv')))


def step_one(datalist, datadir, prec=0.1, angle_prec=5.0, start=-1):
    """ Function that reads through inputdir and uses pymatgen to determine the symmetry group of undeformed structure
    and all it's deformations and writes result to a new .csv file.
    takes symmetry tolerance, angle tolerance and  starting position in datalist as argument (to continue in case of interruptions)
    Defaults are equal to pymatgen default values, and starting from begining of datalist"""

    df = pd.read_csv(datalist, index_col=0, sep=',')
    warnings.filterwarnings("ignore")
    with tqdm.tqdm(total=len(df.index)) as pbar:  # A wrapper that creates nice progress bar
        pbar.set_description("Processing datalist")
        for item in df.index.tolist():
            pbar.update(1)
            for deformation in os.listdir(datadir + str(item)):
                try:
                    poscar = Poscar.from_file(datadir + str(item) + '/' + deformation + '/' + 'POSCAR')
                    structure = poscar.structure
                    df.loc[item, str(deformation)] = str(structure.get_space_group_info(symprec=prec, angle_tolerance=angle_prec))
                except:
                    print("Something is wrong with ", item, ' ', deformation,
                          " and I don't really know why, maybe structure file is not recognized by pymatgen")
                    df.loc[item, str(deformation)] = "pymatgen failure"
    df.to_csv((datalist.replace('.csv', '_sym_check.csv')))


def step_two(datalist, calc_dict):
    """Reads file from previous step and writes TRUE if symmetry of deformation is same
     as symmetry of undeformed, false otherwise"""

    sub_calculations_list = calc_list(calc_dict)
    print(sub_calculations_list)

    df = pd.read_csv(datalist, index_col=0, sep=',')
    warnings.filterwarnings("ignore")
    with tqdm.tqdm(total=len(df.index)) as pbar:  # A wrapper that creates nice progress bar
        pbar.set_description("Processing datalist")
        for item in df.index.tolist():
            pbar.update(1)
            if df.loc[item, 'undeformed'] == "pymatgen failure":
                for sub_calc in sub_calculations_list:
                    df.loc[item, sub_calc] = 'False'
            else:
                for sub_calc in sub_calculations_list:
                    if sub_calc == 'undeformed':
                        pass
                    elif (df.loc[item, sub_calc] == df.loc[item, 'undeformed']) or (pd.isnull(df.loc[item, sub_calc])):
                        df.loc[item, sub_calc] = 'True'
                    else:
                        df.loc[item, sub_calc] = 'False'
                df.loc[item, 'undeformed'] = 'True'

    df.to_csv((datalist.replace('.csv', '_marked.csv')))


def step_three(datalist, calc_dict):
    """Reads file from previous step and drops all entries that have NO problems
    i.e. all deformation symmetries are same as undeformed"""

    sub_calculations_list = calc_list(calc_dict)

    df = pd.read_csv(datalist, index_col=0, sep=',')
    with tqdm.tqdm(total=len(df.index)) as pbar:  # A wrapper that creates nice progress bar
        pbar.set_description("Processing datalist")
        for item in df.index.tolist():
            pbar.update(1)
            result_list = []
            for sub_calc in sub_calculations_list:
                result_list.append(str(df.loc[item, sub_calc]))
            if all(ele == 'True' for ele in result_list):
                df = df.drop([item], axis=0)

    df.to_csv((datalist.replace('_marked.csv', '_warnings.csv')))


def single_test(item):
    """Prints list of symmetries for all possible deformations for single chosen compound. Used for testing."""
    warnings.filterwarnings("ignore")
    for deformation in os.listdir(wdatadir + str(item)):
        try:
            poscar = Poscar.from_file(wdatadir + str(item) + '/' + deformation + '/' + 'POSCAR')
            structure = poscar.structure
            print(deformation, structure.get_space_group_info())
            print(structure.lattice) # uncomment to print lattice matrix
        except:
            print("Something is wrong with ", item, ' ', deformation, " and I don't really know why, maybe structure file is not recognized by pymatgen")


# ------------------------------------------------------------------------------------------------------- #
#  _____                                           _       _____ _             _     _   _                #
# /  __ \                                         | |     /  ___| |           | |   | | | |               #
# | /  \/ ___  _ __ ___  _ __ ___   __ _ _ __   __| |___  \ `--.| |_ __ _ _ __| |_  | |_| | ___ _ __ ___  #
# | |    / _ \| '_ ` _ \| '_ ` _ \ / _` | '_ \ / _` / __|  `--. \ __/ _` | '__| __| |  _  |/ _ \ '__/ _ \ #
# | \__/\ (_) | | | | | | | | | | | (_| | | | | (_| \__ \ /\__/ / || (_| | |  | |_  | | | |  __/ | |  __/ #
#  \____/\___/|_| |_| |_|_| |_| |_|\__,_|_| |_|\__,_|___/ \____/ \__\__,_|_|   \__| \_| |_/\___|_|  \___| #
#                                                                                                         #
# ------------------------------------------------------------------------------------------------------- #
# calculations = {'undeformed': '', 'Applied_Field': '', 'uniaxial': ['0.9', '0.95', '1.05', '1.1'], 'volumetric': ['0.95', '1.05']}

calculations = {'undeformed': '', 'Applied_Field': '', 'uniaxial': ['0.95', '1.05'], 'volumetric': ['0.95', '1.05']}

### Setting which database we work with ###

# # MP
# wdatalist = '../Database/MP/datalist_lattfix_updated_sieved.mag.field_sieved.mag.sites_no.duplicates.csv'
# wdatadir = 'D:/MCES/MP/inputdir/'

# COD
# wdatalist = 'D:/MCES/COD/step0.csv'
# wdatadir = 'D:/MCES/COD/inputdir_step0/'

# # ICSD
wdatalist = 'D:/MCES/ICSD/step0.csv'
wdatadir = 'D:/MCES/ICSD/inputdir_step0/'

# Performing symmetry check:
# Recognising symmetries for this structure and deformations
# step_one(wdatalist, wdatadir)
# Checking which symmetries match with the undeformed structure
# step_two(wdatalist.replace('.csv', '_sym_check.csv'), calculations)
# Creating a datalist containing only entries with some kind of symmetry mismatch
# step_three(wdatalist.replace('.csv', '_sym_check_marked.csv'), calculations)
# Optionally comment following two lines to keep files created during step and step two:
# os.remove(wdatalist.replace('.csv', '_sym_check.csv'))
# os.remove(wdatalist.replace('.csv', '_sym_check_marked.csv'))
#
#single_test('mp-10118')

deformation_check(wdatalist, wdatadir)