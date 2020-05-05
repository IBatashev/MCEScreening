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


def step_one(datalist, prec=0.01, angle_prec=5.0, start=-1):
    """ Function that reads through inputdir and uses pymatgen to determine the symmetry group of undeformed structure
    and all it's deformations and writes result to a new .csv file.
    takes symmetry tolerance, angle tolerance and  starting position in datalist as argument (to continue in case of interruptions)
    Defaults are equal to pymatgen default values, and starting from begining of datalist"""

    df = pd.read_csv(datalist, index_col=0, sep=',')
    warnings.filterwarnings("ignore")
    with tqdm.tqdm(total=len(df.index)) as pbar:  # A wrapper that creates nice progress bar
        pbar.set_description("Processing datalist")
        for item in df.index.tolist():
            if int(item) < start:
                pbar.update(1)
                continue
            else:
                pbar.update(1)
                for deformation in os.listdir(wdatadir + str(item)):
                    try:
                        poscar = Poscar.from_file(wdatadir + str(item) + '/' + deformation + '/' + 'POSCAR')
                        structure = poscar.structure
                        df.loc[item, str(deformation)] = str(structure.get_space_group_info(symprec=prec, angle_tolerance=angle_prec))
                        df.to_csv((datalist.replace('.csv', '_sym_check.csv')))
                    except:
                        print("Something is wrong with ", item, ' ', deformation,
                              " and I don't really know why, maybe structure file is not recognized by pymatgen")
                        df.loc[item, str(deformation)] = "pymatgen failure"


def step_two(datalist):
    """Reads file from previous step and writes TRUE if symmetry of deformation is same
     as symmetry of undeformed, false otherwise"""

    df = pd.read_csv(datalist, index_col=0, sep=',')
    warnings.filterwarnings("ignore")
    with tqdm.tqdm(total=len(df.index)) as pbar:  # A wrapper that creates nice progress bar
        pbar.set_description("Processing datalist")
        for item in df.index.tolist():
            pbar.update(1)
            if df.loc[item, 'undeformed'] == "pymatgen failure":
                df.loc[item, 'a_dec'] = 'False'
                df.loc[item, 'a_inc'] = 'False'
                df.loc[item, 'b_dec'] = 'False'
                df.loc[item, 'b_inc'] = 'False'
                df.loc[item, 'c_dec'] = 'False'
                df.loc[item, 'c_inc'] = 'False'
            else:
                if (df.loc[item, 'a_dec'] == df.loc[item, 'undeformed']) or (pd.isnull(df.loc[item, 'a_dec'])):
                    df.loc[item, 'a_dec'] = 'True'
                else:
                    df.loc[item, 'a_dec'] = 'False'

                if (df.loc[item, 'a_inc'] == df.loc[item, 'undeformed']) or (pd.isnull(df.loc[item, 'a_inc'])):
                    df.loc[item, 'a_inc'] = 'True'
                else:
                    df.loc[item, 'a_inc'] = 'False'

                if (df.loc[item, 'b_dec'] == df.loc[item, 'undeformed']) or (pd.isnull(df.loc[item, 'b_dec'])):
                    df.loc[item, 'b_dec'] = 'True'
                else:
                    df.loc[item, 'b_dec'] = 'False'

                if (df.loc[item, 'b_inc'] == df.loc[item, 'undeformed']) or (pd.isnull(df.loc[item, 'b_inc'])):
                    df.loc[item, 'b_inc'] = 'True'
                else:
                    df.loc[item, 'b_inc'] = 'False'

                if (df.loc[item, 'c_dec'] == df.loc[item, 'undeformed']) or (pd.isnull(df.loc[item, 'c_dec'])):
                    df.loc[item, 'c_dec'] = 'True'
                else:
                    df.loc[item, 'c_dec'] = 'False'

                if (df.loc[item, 'c_inc'] == df.loc[item, 'undeformed']) or (pd.isnull(df.loc[item, 'c_inc'])):
                    df.loc[item, 'c_inc'] = 'True'
                else:
                    df.loc[item, 'c_inc'] = 'False'

    df.to_csv((datalist.replace('.csv', '_marked.csv')))


def step_three(datalist):
    """Reads file from previous step and drops all entries that have NO problems
    i.e. all deformation symmetries are same as undeformed"""

    df = pd.read_csv(datalist, index_col=0, sep=',')
    with tqdm.tqdm(total=len(df.index)) as pbar:  # A wrapper that creates nice progress bar
        pbar.set_description("Processing datalist")
        for item in df.index.tolist():
            pbar.update(1)
            if str(df.loc[item, 'a_dec']) == 'True':
                if str(df.loc[item, 'b_dec']) == 'True':
                    if str(df.loc[item, 'c_dec']) == 'True':
                        if str(df.loc[item, 'a_inc']) == 'True':
                            if str(df.loc[item, 'b_inc']) == 'True':
                                if str(df.loc[item, 'c_inc']) == 'True':
                                    df = df.drop([item], axis=0)
    df.to_csv((datalist.replace('_marked.csv', '_warnings.csv')))


def single_test(item):
    """Prints list of symmetries for all possible deformations for single chosen compound. Used for testing."""
    warnings.filterwarnings("ignore")
    for deformation in os.listdir(wdatadir + str(item)):
        try:
            poscar = Poscar.from_file(wdatadir + str(item) + '/' + deformation + '/' + 'POSCAR')
            structure = poscar.structure

            #print(structure.lattice) # uncomment to print lattice matrix
            print(deformation, structure.get_space_group_info())
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

### Setting which database we work with ###

# wdatalist = '../Database/datalist.csv'
# wdatadir = '../Database/datadir_structure_relaxed/inputdir/'

# wdatalist = '../Database/TESTS/TestDB/datalist_TestDB.csv'
# wdatadir = '../Database/TESTS/TestDB/datadir_TestDB/inputdir/'

wdatalist = '../Database/datalist_updated_sieved.mag.field_sieved.mag.sites_no.duplicates.csv'
wdatadir = 'D:/MCES/inputdir/'

#
# Performing symmetry check:
# Recognising symmetries for this structure and deformations
step_one(wdatalist)
# Checking which symmetries match with the undeformed structure
step_two(wdatalist.replace('.csv', '_sym_check.csv'))
# Creating a datalist containing only entries with some kind of symmetry mismatch
step_three(wdatalist.replace('.csv', '_sym_check_marked.csv'))
# Optionally uncomment following two lines to keep files created during step and step two:
os.remove(wdatalist.replace('.csv', '_sym_check.csv'))
os.remove(wdatalist.replace('.csv', '_sym_check_marked.csv'))

