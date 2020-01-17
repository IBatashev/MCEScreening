import os
import numpy as np

from ase import geometry
import argparse
import shutil
import pandas as pd
import pymatgen
from pymatgen.transformations import standard_transformations  as transform
from pymatgen.io.vasp import Poscar

import INCAR_maker
import POSCAR_maker
import POSCAR_reader
import POTCAR_maker
from tools import POSCAR_reader

### Setting which database we work with ###
wdatadir_structure = '../Database/datadir_structure_relaxed/'
wdatadir_aflow = '../Database/datadir_aflow/'
wdatadir_incars = '../Database/datadir_incars/'
wdatalist = '../Database/datalist.csv'
wdatalist_u = '../Database/datalist_updated_sieved.mag.field_sieved.mag.sites_no.duplicates.csv'


def single_run(def_type, deformation_list):
    """Creates files for a single VASP run"""

    path_to_calc = path + def_type + "/"
    os.mkdir(path_to_calc)
    POSCAR_maker.writer(path_to_calc, poscar_content, deformation_list)
    # INCAR_maker.writer(path_to_calc, poscar_content)
    # POTCAR_maker.writer(path_to_calc, poscar_content)
    # job maker

def undeformed_lattice():
    """Creates all necessary files to run a VASP job for initial undeformed structure"""

    deformation_type = "undeformed"
    deformation = [0,0,0]
    single_run(deformation_type, deformation)

def deformed_lattice(lattice_type, n):
    """Takes lattice type as argument and depending on symmetry creates all necessary files to run VASP jobs for all possible independent deformations.

    For questions concerning this section one should consult bravais_lattice.pdf that contains info on lattice vectors
    for primitive cells of all 14 types. The deformation matrix is determined in a way to certainly add deformation to
    either 'a', 'b' or 'c' parameter and extra unneeded multipliers 'n' would be removed upon multiplication with actual
    lattice matrix in POSCAR_maker.py. e.g.:
     deformation matrix          primitive lattice matrix                                       deformed lattice matrix
    [[n,n,n],[1,1,1],[1,1,1]] * [[3,0,3],[2,0,2],[0,0,1]] = [[3*n,0*n,3*n],[2,0,2],[0,0,1]] = [[3n,0,3n],[2,0,2],[0,0,1]]
    """

    # "CUB"     1 Simple cubic                  OK
    # "FCC"     2 Face centered cubic           OK Some Cubic may be treated wrongly sometimes as tetragonal/orthorhombic
    # "BCC"     3 Body centered cubic           OK
    # "HEX"     4 Hexagonal close packed        mostly good, but 2459 for some reason after streaching a&b is determined as ORCC by ase - what will VASP say? also for 1772 perhaps?
    # "TET"     5 Simple tetragonal             becomes orthorhombic quite often
    # "BCT"     6 Body centered tetragonal      same as above
    # "RHL"     7 Rhombohedral
    # "ORC"     8 Simple orthogonal             #
    # "ORCC"    9 Base centered orthorhombic    # Mostly ok but in some cases becomes MCL
    # "ORCI"    10 Body centered orthorhombic   #
    # "ORCF"    11 Face centered orthorhombic   #
    # "MCL"     12 Simple monoclinic
    # "MCLC"    13 Base centered monoclinic
    # "TRI"     14 Triclinic

    #### OLD WAY I TRIED - NOT GENERAL - works for some, fails for other

    # if lattice_type == "CUB" or lattice_type == "FCC" or lattice_type == "BCC":
    #     deformation_types = 'a'
    #     for deformation_type in deformation_types:
    #         if deformation_type == "a":
    #             deformation_matrix = np.array([[n,n,n],[n,n,n],[n,n,n]])
    #             single_run(deformation_type, deformation_matrix)
    #
    # elif lattice_type == "HEX":
    #     deformation_types = 'ac'
    #     for deformation_type in deformation_types:
    #         if deformation_type == "a":
    #             deformation_matrix = np.array([[n,n,0],[n,n,0],[0,0,1]])
    #             single_run(deformation_type, deformation_matrix)
    #         if deformation_type == "c":
    #             deformation_matrix = np.array([[1,1,0],[1,1,0],[0,0,n]])
    #             single_run(deformation_type, deformation_matrix)
    #
    # elif lattice_type == "RHL":
    #     # Stil confused by how aflow describes rhombohedral lattice BUT for RHL we always have a=b=c so we just need to
    #     # steach everything once (see cubic) and we should be good.
    #     deformation_types = 'a'
    #     for deformation_type in deformation_types:
    #         if deformation_type == "a":
    #             deformation_matrix = np.array([[n,n,n],[n,n,n],[n,n,n]])
    #             single_run(deformation_type, deformation_matrix)
    #
    # elif lattice_type == "TET" or lattice_type == "BCT":
    #     deformation_types = 'ac'
    #     for deformation_type in deformation_types:
    #         if deformation_type == "a":
    #             deformation_matrix = np.array([[n,n,1],[n,n,1],[n,n,1]])
    #             single_run(deformation_type, deformation_matrix)
    #         if deformation_type == "c":
    #             deformation_matrix = np.array([[1,1,n],[1,1,n],[1,1,n]])
    #             single_run(deformation_type, deformation_matrix)
    #
    # elif lattice_type == "ORC" or lattice_type == "ORCC" or lattice_type == "ORCI" or lattice_type == "ORCF":
    #     # the deformation matrixies I now use should properly cover all 4 types of primitive orthorhombic lattices
    #     deformation_types = 'abc'
    #     for deformation_type in deformation_types:
    #         if deformation_type == "a":
    #             deformation_matrix = np.array([[n,1,1],[n,1,1],[n,1,1]])
    #             single_run(deformation_type, deformation_matrix)
    #         if deformation_type == "b":
    #             deformation_matrix = np.array([[1,n,1],[1,n,1],[1,n,1]])
    #             single_run(deformation_type, deformation_matrix)
    #         if deformation_type == "c":
    #             deformation_matrix = np.array([[1,1,n],[1,1,n],[1,1,n]])
    #             single_run(deformation_type, deformation_matrix)
    #
    # elif lattice_type == "MCL" or lattice_type == "MCLC":
    #     deformation_types = 'abc'
    #     for deformation_type in deformation_types:
    #         if deformation_type == "a":
    #             deformation_matrix = np.array([[n,1,1],[n,1,1],[1,1,1]])
    #             single_run(deformation_type, deformation_matrix)
    #         if deformation_type == "b":
    #             deformation_matrix = np.array([[1,n,1],[1,n,1],[1,1,1]])
    #             single_run(deformation_type, deformation_matrix)
    #         if deformation_type == "c":
    #             deformation_matrix = np.array([[1,1,1],[1,1,1],[n,n,n]]) # Aflow for some reason is different for a3 = ccosβx^+csinβz so 'a1' and 'a2' are switched leading to 'a' and 'b' being switched
    #             single_run(deformation_type, deformation_matrix)
    #
    # elif lattice_type == "TRI":
    #     deformation_types = 'abc'
    #     for deformation_type in deformation_types:
    #         if deformation_type == "a":
    #             deformation_matrix = np.array([[n,0,0],[1,1,0],[1,1,1]])
    #             single_run(deformation_type, deformation_matrix)
    #         if deformation_type == "b":
    #             deformation_matrix = np.array([[1,0,0],[n,n,0],[1,1,1]])
    #             single_run(deformation_type, deformation_matrix)
    #         if deformation_type == "c":
    #             deformation_matrix = np.array([[1,0,0],[1,1,0],[n,n,n]])
    #             single_run(deformation_type, deformation_matrix)


    if lattice_type == "CUB" or lattice_type == "FCC" or lattice_type == "BCC":
        deformation_types = 'a'
        for deformation_type in deformation_types:
            if deformation_type == "a":
                deformation = [n, n, n]
                single_run(deformation_type, deformation)

    elif lattice_type == "HEX":
        deformation_types = 'ac'
        for deformation_type in deformation_types:
            if deformation_type == "a":
                deformation = [n, n, 0]
                single_run(deformation_type, deformation)
            if deformation_type == "c":
                deformation = [0, 0, n]
                single_run(deformation_type, deformation)

    elif lattice_type == "RHL":
        deformation_types = 'a'
        for deformation_type in deformation_types:
            if deformation_type == "a":
                deformation = [n, n, n]
                single_run(deformation_type, deformation)

    elif lattice_type == "TET" or lattice_type == "BCT":
        deformation_types = 'ac'
        for deformation_type in deformation_types:
            if deformation_type == "a":
                deformation = [n, n, 0]
                single_run(deformation_type, deformation)
            if deformation_type == "c":
                deformation = [0, 0, n]
                single_run(deformation_type, deformation)

    elif lattice_type == "ORC" or lattice_type == "ORCC" or lattice_type == "ORCI" or lattice_type == "ORCF":
        deformation_types = 'abc'
        for deformation_type in deformation_types:
            if deformation_type == "a":
                deformation = [n, 0, 0]
                single_run(deformation_type, deformation)
            if deformation_type == "b":
                deformation = [0, n, 0]
                single_run(deformation_type, deformation)
            if deformation_type == "c":
                deformation = [0, 0, n]
                single_run(deformation_type, deformation)

    elif lattice_type == "MCL" or lattice_type == "MCLC":
        deformation_types = 'abc'
        for deformation_type in deformation_types:
            if deformation_type == "a":
                deformation = [n, 0, 0]
                single_run(deformation_type, deformation)
            if deformation_type == "b":
                deformation = [0, n, 0]
                single_run(deformation_type, deformation)
            if deformation_type == "c":
                deformation = [0, 0, n]
                single_run(deformation_type, deformation)

    elif lattice_type == "TRI":
        deformation_types = 'abc'
        for deformation_type in deformation_types:
            if deformation_type == "a":
                deformation = [n, 0, 0]
                single_run(deformation_type, deformation)
            if deformation_type == "b":
                deformation = [0, n, 0]
                single_run(deformation_type, deformation)
            if deformation_type == "c":
                deformation = [0, 0, n]
                single_run(deformation_type, deformation)

    else: print("Error! Unknown Lattice type for " + str(ID))

def parseArguments():
    """Function for parsing input argumens necessary for creator.py to work.
    We expect to at least get ID (name of POSCAR) to work with.
    By default deformation coefficient is set to 1.2"""

    parser = argparse.ArgumentParser()    # Create argument parser
    # Positional mandatory arguments
    parser.add_argument("ID", help="Name of POSCAR file", type=str)
    # Optional arguments
    parser.add_argument("-d", "--deformation", help="Deformation Coefficient", type=float, default=1.2)
    # Parse arguments
    args = parser.parse_args()
    return args


#--------------------- MAIN PART STARTS HERE ---------------------#
# args = parseArguments() # Get input arguments:
# ID = args.ID            # ID of the structure we are going to work with
# n = args.deformation    # deformation coefficient - an optional input argument
# path = '../'+ID+'/'     # Prepare folder folder where all subfolders for a single ID
# os.makedirs(path)       # will be created/executed/cleaned - working directory

ID = 37                  #
def_factor  = 0.2                     #
path = '../Sample_checks/'+str(ID)+'/'     # Prepare folder folder where all subfolders for a single ID
if os.path.exists(path):    # TEMP for local tests
    shutil.rmtree(path)     #
os.makedirs(path)           #

### Take all information from poscar as list of strings:
structure_file = wdatadir_structure + str(ID)
poscar_content = POSCAR_reader.read(structure_file)
df = pd.read_csv(wdatalist_u, index_col=0, sep=',')
latt_type = (df.loc[ID, 'lattice_system'])                     # not sure what is more reliable or even if there is any difference between bravais type or just name of system


### Use ase to get Bravais lattice type from lattice given in POSCAR - ase seems to be a bit more adequate but thismatch is only 27 out of ~1200
lattice_matrix = np.zeros([3, 3])
for i in range(2, 5):
    lattice_matrix[i - 2, :] = np.fromstring(poscar_content[i], dtype=np.float, sep=' ')
read_lattice_matrix = geometry.Cell.new(lattice_matrix)
lat_type = str(read_lattice_matrix.get_bravais_lattice()).split("(")[0]

if def_factor == 0:            # Check if we want deformation to happen
    undeformed_lattice()       # Create a folder for a job without any deformation and run it
else:
    undeformed_lattice()
    deformed_lattice(lat_type, def_factor, )  # Create the proper number of folders for all possible deformations according to Bravais lattice type and run them









 # kk = geometry.cell_to_cellpar(read_lattice_matrix)
    # struct = geometry.cellpar_to_cell(kk)
    # struct[0] = 1*struct[0]
    # struct[1] = 1 * struct[1]
    # struct[2] = 1 * struct[2]
    #
    # kkk = geometry.cell_to_cellpar(struct)
    # return lattice_matrix, kk, struct, kkk

# inp = '../Database/datalist_updated_sieved.mag.field_sieved.mag.sites_no.duplicates_tests.bravais.csv'
# df = pd.read_csv(inp, index_col=0, sep=',')
#
# for item in df.index.tolist():
#     # ID = item
#     # df.loc[item, 'ase_bravais'] = creator(item, 1.0)
#     # print(df.loc[item, 'ase_bravais'], df.loc[item, 'bravais_lattice'])
#     if (df.loc[item, 'ase_bravais'] == df.loc[item, 'bravais_lattice']):
#         df = df.drop([item])
# df.to_csv(inp.replace(".csv", '_tests.bravais.dropped.' + '.csv'))

# IDN = 1724
#
# print(creator(IDN)[0], '\n')
# print(creator(IDN)[1], '\n')
# print(creator(IDN)[2], '\n')
# print(creator(IDN)[3])


# poscar = Poscar.from_file("POSCAR")
def deformed_tester(ID, n):
    df = pd.read_csv(wdatalist_u, index_col=0, sep=',')
    lattice_type = (df.loc[ID, 'lattice_system'])
    poscar = POSCAR_reader.read(wdatadir_structure + str(ID))
    poscar_string = ''.join(poscar)
    poscar = Poscar.from_string(poscar_string)
    structure = poscar.structure
    print(structure.lattice)
    print(structure)
    #
    print('after')


    if lattice_type == 'cubic':
        deformation = [n, n, n]
        deformation_type ='a'
        single_run(deformation_type, deformation)


    structure.apply_strain((n, 0, 0))
    print(structure.lattice)
    print(structure)

# deformed_tester(37, 0.2)
