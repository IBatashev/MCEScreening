import pandas as pd
from pymatgen.io.vasp import Poscar
import os
import tqdm
import warnings
wdatalist = '../Database/datalist.csv'
wdatadir = '../Database/datadir_structure_relaxed/inputdir/'



def one():
    df = pd.read_csv(wdatalist, index_col=0, sep=',')
    warnings.filterwarnings("ignore")
    with tqdm.tqdm(total=len(df.index)) as pbar:  # A wrapper that creates nice progress bar
        pbar.set_description("Processing datalist")
        for item in df.index.tolist():
            if int(item) <= 7513:
                pbar.update(1)
                continue
            else:
                pbar.update(1)


                for deformation in os.listdir(wdatadir + str(item)):
                    poscar = Poscar.from_file(wdatadir + str(item) + '/' + deformation + '/' + 'POSCAR')
                    structure = poscar.structure
                    # print(structure)
                    df.loc[item, str(deformation)] = str(structure.get_space_group_info())
                df.to_csv((wdatalist.replace('.csv', '_sym_check.csv')))


def two():
    """Removes entries from datalist according to selected sieving criteria"""

    df = pd.read_csv('../Database/datalist_sym_check_marked.csv', index_col=0, sep=',')
    with tqdm.tqdm(total=len(df.index)) as pbar:  # A wrapper that creates nice progress bar
        pbar.set_description("Processing datalist")
        for item in df.index.tolist():
            pbar.update(1)
            if str(df.loc[item, 'a_check']) == 'True' :
                if str(df.loc[item, 'b_check'])== 'True' :
                    if str(df.loc[item, 'c_check']) == 'True' :
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


# deformed_tester(7517, 0.2, 'b')
two()