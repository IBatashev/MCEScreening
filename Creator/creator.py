import os
import numpy as np

from ase import geometry
import argparse
import shutil
import pandas as pd
import tqdm
import pymatgen
from pymatgen.transformations import standard_transformations  as transform
from pymatgen.io.vasp import Poscar
import warnings
import INCAR_maker
import POSCAR_maker
import POSCAR_reader
import POTCAR_maker
# from tools import POSCAR_reader

### Setting which database we work with ###
wdatadir_structure = '../Database/datadir_structure_relaxed/'
wdatadir_aflow = '../Database/datadir_aflow/'
wdatadir_incars = '../Database/datadir_incars/'
wdatalist = '../Database/datalist.csv'

wdatalist = '../Database/TestDB/datalist_TestDB.csv'
wdatadir_structure = '../Database/TestDB/datadir_TestDB/'

def single_run(def_type, deformation_list):
    """Creates files for a single VASP run"""

    path_to_calc = path + def_type + "/"
    os.mkdir(path_to_calc)
    POSCAR_maker.writer(path_to_calc, poscar_content, deformation_list)
    INCAR_maker.writer(path_to_calc, poscar_content, elem_list)
    POTCAR_maker.writer(path_to_calc, elem_list)


def undeformed_lattice():
    """Creates all necessary files to run a VASP job for initial undeformed structure"""

    deformation_type = "undeformed"
    deformation = [[1,1,1],[1,1,1],[1,1,1]]
    single_run(deformation_type, deformation)


def deformed_lattice(lattice_type, n):
    """Takes lattice type as argument and depending on symmetry creates all necessary files to run VASP jobs for all
    possible independent deformations.

    Depending on symmetry create deformation list [n, n, n ] which determines weather a, b, c are deformed. This list is
    later passed to POSCAR_maker where pymatgen is used for adjusting the lattice matrix.
    E.g. deformation [0.2, 0, 0] means that we have a parameter increased by 20% and no changes to b and c.
    Use negative numbers for compression.
    """

   # # "CUB"     1 Simple cubic                  OK
   # # "FCC"     2 Face centered cubic           OK
   # # "BCC"     3 Body centered cubic           OK
   # # "HEX"     4 Hexagonal close packed        Mostly OK, vefy few are identified wrongly by pymatgen, fixed by increasing symprec  - test with vasp
   # # "TET"     5 Simple tetragonal             seems ok - only 3/8603 are wrongly identified by pymatgen may test them with vasp
   # # "BCT"     6 Body centered tetragonal      shit, fix needed - problem with pymatgen deformation method, need my own.
   # # "RHL"     7 Rhombohedral                  OK, one spacegroup missinterpretation and it still belongs to rhombohedral symmetry
   # # "ORC"     8 Simple orthogonal             OK
   # # "ORCC"    9 Base centered orthorhombic    Shit
   # # "ORCI"    10 Body centered orthorhombic   shit
   # # "ORCF"    11 Face centered orthorhombic   shit
   # # "MCL"     12 Simple monoclinic            few mistakes, need to see why, maby same story as hex?
   # # "MCLC"    13 Base centered monoclinic     shit
   # # "TRI"     14 Triclinic                    shit but looks working just because this is the lowest symmetry - cannot fall worse


    if lattice_type == "CUB" or lattice_type == "FCC" or lattice_type == "BCC":
        deformation_types = 'a'
        for deformation_type in deformation_types:
            if deformation_type == "a":
                deformation_matrix = np.array([[n,n,n],[n,n,n],[n,n,n]])
                single_run(deformation_type, deformation_matrix)

    elif lattice_type == "HEX":
        deformation_types = 'acV'
        for deformation_type in deformation_types:
            if deformation_type == "a":
                deformation_matrix = np.array([[n,n,0],[n,n,0],[0,0,1]])
                single_run(deformation_type, deformation_matrix)
            if deformation_type == "c":
                deformation_matrix = np.array([[1,1,0],[1,1,0],[0,0,n]])
                single_run(deformation_type, deformation_matrix)
            if deformation_type == "V":
                deformation_matrix = np.array([[n,n,n],[n,n,n],[n,n,n]])
                single_run(deformation_type, deformation_matrix)

    elif lattice_type == "RHL":
        # Stil confused by how aflow describes rhombohedral lattice BUT for RHL we always have a=b=c so we just need to
        # steach everything once (see cubic) and we should be good.
        deformation_types = 'a'
        for deformation_type in deformation_types:
            if deformation_type == "a":
                deformation_matrix = np.array([[n,n,n],[n,n,n],[n,n,n]])
                single_run(deformation_type, deformation_matrix)

    elif lattice_type == "TET" or lattice_type == "BCT":
        deformation_types = 'acV'
        for deformation_type in deformation_types:
            if deformation_type == "a":
                deformation_matrix = np.array([[n,n,1],[n,n,1],[n,n,1]])
                single_run(deformation_type, deformation_matrix)
            if deformation_type == "c":
                deformation_matrix = np.array([[1,1,n],[1,1,n],[1,1,n]])
                single_run(deformation_type, deformation_matrix)
            if deformation_type == "V":
                deformation_matrix = np.array([[n,n,n],[n,n,n],[n,n,n]])
                single_run(deformation_type, deformation_matrix)

    elif lattice_type == "ORC" or lattice_type == "ORCC" or lattice_type == "ORCI" or lattice_type == "ORCF":
        # the deformation matrixies I now use should properly cover all 4 types of primitive orthorhombic lattices
        deformation_types = 'abcV'
        for deformation_type in deformation_types:
            if deformation_type == "a":
                deformation_matrix = np.array([[n,1,1],[n,1,1],[n,1,1]])
                single_run(deformation_type, deformation_matrix)
            if deformation_type == "b":
                deformation_matrix = np.array([[1,n,1],[1,n,1],[1,n,1]])
                single_run(deformation_type, deformation_matrix)
            if deformation_type == "c":
                deformation_matrix = np.array([[1,1,n],[1,1,n],[1,1,n]])
                single_run(deformation_type, deformation_matrix)
            if deformation_type == "V":
                deformation_matrix = np.array([[n, n, n], [n, n, n], [n, n, n]])
                single_run(deformation_type, deformation_matrix)

    elif lattice_type == "MCL" or lattice_type == "MCLC":
        deformation_types = 'abcV'
        for deformation_type in deformation_types:
            if deformation_type == "a":
                deformation_matrix = np.array([[n,1,1], [n,1,1], [1,1,1]])
                single_run(deformation_type, deformation_matrix)
            if deformation_type == "b":
                deformation_matrix = np.array([[1,n,1], [1,n,1], [1,1,1]])
                single_run(deformation_type, deformation_matrix)
            if deformation_type == "c":
                deformation_matrix = np.array([[1,1,1], [1,1,1], [n,n,n]]) # Aflow for some reason is different for a3 = ccosβx^+csinβz so 'a1' and 'a2' are switched leading to 'a' and 'b' being switched
                single_run(deformation_type, deformation_matrix)
            if deformation_type == "V":
                deformation_matrix = np.array([[n, n, n], [n, n, n], [n, n, n]])
                single_run(deformation_type, deformation_matrix)

    elif lattice_type == "TRI":
        deformation_types = 'abcV'
        for deformation_type in deformation_types:
            if deformation_type == "a":
                deformation_matrix = np.array([[n,0,0],[1,1,0],[1,1,1]])
                single_run(deformation_type, deformation_matrix)
            if deformation_type == "b":
                deformation_matrix = np.array([[1,0,0],[n,n,0],[1,1,1]])
                single_run(deformation_type, deformation_matrix)
            if deformation_type == "c":
                deformation_matrix = np.array([[1,0,0],[1,1,0],[n,n,n]])
                single_run(deformation_type, deformation_matrix)
            if deformation_type == "V":
                deformation_matrix = np.array([[n,n,n],[n,n,n],[n,n,n]])
                single_run(deformation_type, deformation_matrix)

    # if lattice_type == "CUB" or lattice_type == "FCC" or lattice_type == "BCC":
    #     deformation_types = 'a'
    #     for deformation_type in deformation_types:
    #         if deformation_type == "a":
    #             deformation = [n, n, n]
    #             single_run(deformation_type, deformation)
    #
    # elif lattice_type == "HEX":
    #     deformation_types = 'ac'
    #     for deformation_type in deformation_types:
    #         if deformation_type == "a":
    #             deformation = [n, n, 0]
    #             single_run(deformation_type, deformation)
    #         if deformation_type == "c":
    #             deformation = [0, 0, n]
    #             single_run(deformation_type, deformation)
    #
    # elif lattice_type == "RHL":
    #     deformation_types = 'a'
    #     for deformation_type in deformation_types:
    #         if deformation_type == "a":
    #             deformation = [n, n, n]
    #             single_run(deformation_type, deformation)
    #
    # elif lattice_type == "tetragonal":
    #     deformation_types = 'ac'
    #     for deformation_type in deformation_types:
    #         if deformation_type == "a":
    #             deformation = [n, n, 0]
    #             single_run(deformation_type, deformation)
    #         if deformation_type == "c":
    #             deformation = [0, 0, n]
    #             single_run(deformation_type, deformation)
    #
    # elif lattice_type == "orthorhombic":
    #     deformation_types = 'abc'
    #     for deformation_type in deformation_types:
    #         if deformation_type == "a":
    #             deformation = [n, 0, 0]
    #             single_run(deformation_type, deformation)
    #         if deformation_type == "b":
    #             deformation = [0, n, 0]
    #             single_run(deformation_type, deformation)
    #         if deformation_type == "c":
    #             deformation = [0, 0, n]
    #             single_run(deformation_type, deformation)
    #
    # elif lattice_type == "monoclinic":
    #     deformation_types = 'abc'
    #     for deformation_type in deformation_types:
    #         if deformation_type == "a":
    #             deformation = [n, 0, 0]
    #             single_run(deformation_type, deformation)
    #         if deformation_type == "b":
    #             deformation = [0, n, 0]
    #             single_run(deformation_type, deformation)
    #         if deformation_type == "c":
    #             deformation = [0, 0, n]
    #             single_run(deformation_type, deformation)
    #
    # elif lattice_type == "triclinic":
    #     deformation_types = 'abc'
    #     for deformation_type in deformation_types:
    #         if deformation_type == "a":
    #             deformation = [n, 0, 0]
    #             single_run(deformation_type, deformation)
    #         if deformation_type == "b":
    #             deformation = [0, n, 0]
    #             single_run(deformation_type, deformation)
    #         if deformation_type == "c":
    #             deformation = [0, 0, n]
    #             single_run(deformation_type, deformation)

    else:
        print("Error! Unknown Lattice type for " + str(ID))


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


# --------------------- MAIN PART STARTS HERE ---------------------#
# args = parseArguments() # Get input arguments:
# ID = args.ID            # ID of the structure we are going to work with
# n = args.deformation    # deformation coefficient - an optional input argument
# path = '../'+ID+'/'     # Prepare folder folder where all subfolders for a single ID
# os.makedirs(path)       # will be created/executed/cleaned - working directory

warnings.filterwarnings("ignore") # slightly dirty, but simple way to disable pymatgen complaining about lack
                                  # of element labels in POSCARs from aflowlib "UserWarning: Elements in POSCAR cannot be determined."
                                  # Warning ! This disables _ALL_ warnings.



df = pd.read_csv(wdatalist, index_col=0, sep=',')

with tqdm.tqdm(total=len(df.index)) as pbar:  # A wrapper that creates nice progress bar
    pbar.set_description("Processing datalist")
    for item in df.index.tolist():
        pbar.update(1)
        ID = str(item)
        def_factor = 1.0
        path = wdatadir_structure + 'inputdir/'+str(ID)+'/'     # Prepare folder folder where all subfolders for a single ID
        if os.path.exists(path):    # !WARNING! Overwrites existing path!
            shutil.rmtree(path)     #
        os.makedirs(path)           #

        ### Take all information from poscar as list of strings:
        structure_file = wdatadir_structure + str(ID)
        poscar_content = POSCAR_reader.read(structure_file)
        # df = pd.read_csv(wdatalist_u, index_col=0, sep=',')
        lat_type = (df.loc[item, 'Bravais_lattice'])
        elem_list = (df.loc[item, 'species']).replace(';', ',').replace("'", "").strip('][').split(', ')

        if def_factor == 1:            # Check if we want deformation to happen
            undeformed_lattice()       # Create a folder for a job without any deformation
        else:
            undeformed_lattice()
            deformed_lattice(lat_type, def_factor)  # Create the proper number of folders for all possible deformations according to Bravais lattice type
