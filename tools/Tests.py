import pandas as pd
from pymatgen.io.vasp import Poscar
import os
import tqdm
import warnings

wdatalist = '../Database/datalist.csv'
wdatadir = '../Database/datadir_structure_relaxed/inputdir/'

wdatalist = '../Database/TestDB/datalist_TestDB.csv'
wdatadir = '../Database/TestDB/datadir_TestDB/inputdir/'

# This module should be used to make sure we get symmetries in all our structures properly recognized by vasp
# checks are perforemed in several steps

def one(prec=0.01, angle_prec=5.0, start=-1):
    """ Function that reads through inputdir and uses pymatgen to determine the symmetry group of undeformed structure
    and all it's deformations and writes result to a new .csv file.
    takes symmetry tolerance, angle tolerance and  starting position in datalist as argument (to continue in case of interruptions)
    Defaults are equal to pymatgen default values, and starting from begining of datalist"""

    df = pd.read_csv(wdatalist, index_col=0, sep=',')
    warnings.filterwarnings("ignore")
    with tqdm.tqdm(total=len(df.index)) as pbar:  # A wrapper that creates nice progress bar
        pbar.set_description("Processing datalist")
        for item in df.index.tolist():
            if int(item) < start:
                pbar.update(1)          # to skip already done part, note sometimes seem to behave a bit strange perhaps because ID is used instead of counter
                continue
            else:
                pbar.update(1)
                for deformation in os.listdir(wdatadir + str(item)):
                    poscar = Poscar.from_file(wdatadir + str(item) + '/' + deformation + '/' + 'POSCAR')
                    structure = poscar.structure
                    df.loc[item, str(deformation)] = str(structure.get_space_group_info(symprec=prec, angle_tolerance=angle_prec))
                    df.to_csv((wdatalist.replace('.csv', '_sym_check.csv')))


def two():
    """read file from previous step and writes TRUE if symmetry of deformation is same as symmetry of undeformed, false otherwise"""

    # for now is done in excel should add this function later or just combine whith step three

def three():
    """Drops all entries that have NO problems i.e. all deformation symmetries are same as undeformed"""

    df = pd.read_csv('../Database/datalist_sym_check_marked.csv', index_col=0, sep=',')
    with tqdm.tqdm(total=len(df.index)) as pbar:  # A wrapper that creates nice progress bar
        pbar.set_description("Processing datalist")
        for item in df.index.tolist():
            pbar.update(1)
            if str(df.loc[item, 'a_check']) == 'True' :
                if str(df.loc[item, 'b_check'])== 'True' :
                    if str(df.loc[item, 'c_check']) == 'True' :
                        if str(df.loc[item, 'V_check']) == 'True' :
                            df = df.drop([item], axis=0)
    df.to_csv('../Database/datalist_sym_check_marked_sieved.csv')

def deformed_tester(ID, n, deformation):
    df = pd.read_csv(wdatalist, index_col=0, sep=',')
    # lattice_type = (df.loc[ID, 'lattice_system'])
    # poscar = POSCAR_reader.read(wdatadir + str(ID))
    poscar1 = Poscar.from_file(wdatadir + str(ID) + '/' + 'undeformed' + '/' + "POSCAR")
    # poscar_string = ''.join(poscar)
    # poscar = Poscar.from_string(poscar_string)
    structure1 = poscar1.structure

    # poscar2 = POSCAR_reader.read('../Creator/5565')
    # poscar_string2 = ''.join(poscar2)
    poscar2 = Poscar.from_file(wdatadir + str(ID) + '/' + deformation + '/' + "POSCAR")
    structure2 = poscar2.structure

    print(structure1.lattice)
    print(structure1)
    print(structure1.get_space_group_info())
    #
    print('after')
    # if lattice_type == 'cubic':
    #     deformation = [n, n, n]
    #     deformation_type ='a'
    #     single_run(deformation_type, deformation)
    #
    # structure1.apply_strain((0, 0, 0))
    print(structure2.lattice)
    print(structure2)
    print(structure2.get_space_group_info(symprec=0.01))


def one_single(item):
    df = pd.read_csv(wdatalist, index_col=0, sep=',')
    warnings.filterwarnings("ignore")
    for deformation in os.listdir(wdatadir + str(item)):
        poscar = Poscar.from_file(wdatadir + str(item) + '/' + deformation + '/' + 'POSCAR')
        # print(poscar)

        structure = poscar.structure
        #print(structure.lattice)
        print(deformation, structure.get_space_group_info())
        # df.loc[item, str(deformation)] = str(structure.get_space_group_info())


# deformed_tester(7517, 0.2, 'b')
one()
# one_single(2692)