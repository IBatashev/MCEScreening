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
    deformation = [[1, 1, 1], [1, 1, 1], [1, 1, 1]]
    single_run(deformation_type, deformation)


def deformed_lattice(lattice_type, n):
    """Takes lattice type as argument and depending on symmetry creates all necessary files to run VASP jobs for all
    possible independent deformations.
    Depending on symmetry create a 3x3 deformation matrix which determines weather a, b, c are deformed. This matrix is
    later passed to POSCAR_maker where it is used for adjusting the lattice matrix.
    Deformations are done both for increase and decrease of each lattice parameter.
    """

    # "CUB"     1 Simple cubic
    # "FCC"     2 Face centered cubic
    # "BCC"     3 Body centered cubic
    # "HEX"     4 Hexagonal close packed
    # "TET"     5 Simple tetragonal
    # "BCT"     6 Body centered tetragonal
    # "RHL"     7 Rhombohedral
    # "ORC"     8 Simple orthogonal
    # "ORCC"    9 Base centered orthorhombic
    # "ORCI"    10 Body centered orthorhombic
    # "ORCF"    11 Face centered orthorhombic
    # "MCL"     12 Simple monoclinic
    # "MCLC"    13 Base centered monoclinic
    # "TRI"     14 Triclinic

    if lattice_type == "CUB" or lattice_type == "FCC" or lattice_type == "BCC":
        deformation_types = ['a_dec', 'a_inc']
        for deformation_type in deformation_types:
            if deformation_type == "a_dec":
                deformation_matrix = np.array([[1-n, 1-n, 1-n], [1-n, 1-n, 1-n], [1-n, 1-n, 1-n]])
                single_run(deformation_type, deformation_matrix)
            if deformation_type == "a_inc":
                deformation_matrix = np.array([[1+n, 1+n, 1+n], [1+n, 1+n, 1+n], [1+n, 1+n, 1+n]])
                single_run(deformation_type, deformation_matrix)

    elif lattice_type == "HEX":
        deformation_types = ['a_dec', 'c_dec', 'a_inc', 'c_inc']
        for deformation_type in deformation_types:
            if deformation_type == "a_dec":
                deformation_matrix = np.array([[1-n, 1-n, 0], [1-n, 1-n, 0], [0, 0, 1]])
                single_run(deformation_type, deformation_matrix)
            if deformation_type == "c_dec":
                deformation_matrix = np.array([[1, 1, 0], [1, 1, 0], [0, 0, 1-n]])
                single_run(deformation_type, deformation_matrix)
            if deformation_type == "a_inc":
                deformation_matrix = np.array([[1+n, 1+n, 0], [1+n, 1+n, 0], [0, 0, 1]])
                single_run(deformation_type, deformation_matrix)
            if deformation_type == "c_inc":
                deformation_matrix = np.array([[1, 1, 0], [1, 1, 0], [0, 0, 1+n]])
                single_run(deformation_type, deformation_matrix)

    elif lattice_type == "RHL":
        # Stil confused by how aflow describes rhombohedral lattice BUT for RHL we always have a=b=c so we just need to
        # steach everything once (see cubic) and we should be good.
        deformation_types = ['a_dec', 'a_inc']
        for deformation_type in deformation_types:
            if deformation_type == "a_dec":
                deformation_matrix = np.array([[1-n, 1-n, 1-n], [1-n, 1-n, 1-n], [1-n, 1-n, 1-n]])
                single_run(deformation_type, deformation_matrix)
            if deformation_type == "a_inc":
                deformation_matrix = np.array([[1+n, 1+n, 1+n], [1+n, 1+n, 1+n], [1+n, 1+n, 1+n]])
                single_run(deformation_type, deformation_matrix)

    elif lattice_type == "TET" or lattice_type == "BCT":
        deformation_types = ['a_dec', 'c_dec', 'a_inc', 'c_inc']
        for deformation_type in deformation_types:
            if deformation_type == "a_dec":
                deformation_matrix = np.array([[1-n, 1-n, 1], [1-n, 1-n, 1], [1-n, 1-n, 1]])
                single_run(deformation_type, deformation_matrix)
            if deformation_type == "c_dec":
                deformation_matrix = np.array([[1, 1, 1-n], [1, 1, 1-n], [1, 1, 1-n]])
                single_run(deformation_type, deformation_matrix)
            if deformation_type == "a_inc":
                deformation_matrix = np.array([[1+n, 1+n, 1], [1+n, 1+n, 1], [1+n, 1+n, 1]])
                single_run(deformation_type, deformation_matrix)
            if deformation_type == "c_inc":
                deformation_matrix = np.array([[1, 1, 1+n], [1, 1, 1+n], [1, 1, 1+n]])
                single_run(deformation_type, deformation_matrix)

    elif lattice_type == "ORC" or lattice_type == "ORCC" or lattice_type == "ORCI" or lattice_type == "ORCF":
        deformation_types = ['a_dec', 'b_dec', 'c_dec', 'a_inc', 'b_inc', 'c_inc']
        for deformation_type in deformation_types:
            if deformation_type == "a_dec":
                deformation_matrix = np.array([[1-n, 1, 1], [1-n, 1, 1], [1-n, 1, 1]])
                single_run(deformation_type, deformation_matrix)
            if deformation_type == "b_dec":
                deformation_matrix = np.array([[1, 1-n, 1], [1, 1-n, 1], [1, 1-n, 1]])
                single_run(deformation_type, deformation_matrix)
            if deformation_type == "c_dec":
                deformation_matrix = np.array([[1, 1, 1-n], [1, 1, 1-n], [1, 1, 1-n]])
                single_run(deformation_type, deformation_matrix)
            if deformation_type == "a_inc":
                deformation_matrix = np.array([[1+n, 1, 1], [1+n, 1, 1], [1+n, 1, 1]])
                single_run(deformation_type, deformation_matrix)
            if deformation_type == "b_inc":
                deformation_matrix = np.array([[1, 1+n, 1], [1, 1+n, 1], [1, 1+n, 1]])
                single_run(deformation_type, deformation_matrix)
            if deformation_type == "c_inc":
                deformation_matrix = np.array([[1, 1, 1+n], [1, 1, 1+n], [1, 1, 1+n]])
                single_run(deformation_type, deformation_matrix)

    elif lattice_type == "MCL" or lattice_type == "MCLC":
        deformation_types = ['a_dec', 'b_dec', 'c_dec', 'a_inc', 'b_inc', 'c_inc']
        for deformation_type in deformation_types:
            if deformation_type == "a_dec":
                deformation_matrix = np.array([[1-n, 1, 1], [1-n, 1, 1], [1, 1, 1]])
                single_run(deformation_type, deformation_matrix)
            if deformation_type == "b_dec":
                deformation_matrix = np.array([[1, 1-n, 1], [1, 1-n, 1], [1, 1, 1]])
                single_run(deformation_type, deformation_matrix)
            if deformation_type == "c_dec":
                deformation_matrix = np.array([[1, 1, 1], [1, 1, 1], [1-n, 1-n, 1-n]]) # Aflow for some reason is different for a3 = ccosβx^+csinβz so 'a1' and 'a2' are switched leading to 'a' and 'b' being switched
                single_run(deformation_type, deformation_matrix)
            if deformation_type == "a_inc":
                deformation_matrix = np.array([[1+n, 1, 1], [1+n, 1, 1], [1, 1, 1]])
                single_run(deformation_type, deformation_matrix)
            if deformation_type == "b_inc":
                deformation_matrix = np.array([[1, 1+n, 1], [1, 1+n, 1], [1, 1, 1]])
                single_run(deformation_type, deformation_matrix)
            if deformation_type == "c_inc":
                deformation_matrix = np.array([[1, 1, 1], [1, 1, 1], [1+n, 1+n, 1+n]]) # Aflow for some reason is different for a3 = ccosβx^+csinβz so 'a1' and 'a2' are switched leading to 'a' and 'b' being switched
                single_run(deformation_type, deformation_matrix)

    elif lattice_type == "TRI":
        deformation_types = ['a_dec', 'b_dec', 'c_dec', 'a_inc', 'b_inc', 'c_inc']
        for deformation_type in deformation_types:
            if deformation_type == "a_dec":
                deformation_matrix = np.array([[1-n, 0, 0], [1, 1, 0], [1, 1, 1]])
                single_run(deformation_type, deformation_matrix)
            if deformation_type == "b_dec":
                deformation_matrix = np.array([[1, 0, 0], [1-n, 1-n, 0], [1, 1, 1]])
                single_run(deformation_type, deformation_matrix)
            if deformation_type == "c_dec":
                deformation_matrix = np.array([[1, 0, 0], [1, 1, 0], [1-n, 1-n, 1-n]])
                single_run(deformation_type, deformation_matrix)
            if deformation_type == "a_inc":
                deformation_matrix = np.array([[1+n, 0, 0], [1, 1, 0], [1, 1, 1]])
                single_run(deformation_type, deformation_matrix)
            if deformation_type == "b_inc":
                deformation_matrix = np.array([[1, 0, 0], [1+n, 1+n, 0], [1, 1, 1]])
                single_run(deformation_type, deformation_matrix)
            if deformation_type == "c_inc":
                deformation_matrix = np.array([[1, 0, 0], [1, 1, 0], [1+n, 1+n, 1+n]])
                single_run(deformation_type, deformation_matrix)

    else:
        print("Error! Unknown Lattice type for " + str(ID))


def parseArguments():
    """Function for parsing input arguments necessary for creator.py to work.
    We expect to at least get ID (name of POSCAR) to work with.
    By default deformation coefficient is set to 5%
    THIS FUNCTION IS CURRENTLY UNUSED AND OUTDATED"""

    parser = argparse.ArgumentParser()    # Create argument parser
    # Positional mandatory arguments
    parser.add_argument("ID", help="Name of POSCAR file", type=str)
    # Optional arguments
    parser.add_argument("-d", "--deformation", help="Deformation Coefficient", type=float, default=0.05)
    # Parse arguments
    args = parser.parse_args()
    return args


# ------------------------------------------------------------------------------------------------------- #
#  _____                                           _       _____ _             _     _   _                #
# /  __ \                                         | |     /  ___| |           | |   | | | |               #
# | /  \/ ___  _ __ ___  _ __ ___   __ _ _ __   __| |___  \ `--.| |_ __ _ _ __| |_  | |_| | ___ _ __ ___  #
# | |    / _ \| '_ ` _ \| '_ ` _ \ / _` | '_ \ / _` / __|  `--. \ __/ _` | '__| __| |  _  |/ _ \ '__/ _ \ #
# | \__/\ (_) | | | | | | | | | | | (_| | | | | (_| \__ \ /\__/ / || (_| | |  | |_  | | | |  __/ | |  __/ #
#  \____/\___/|_| |_| |_|_| |_| |_|\__,_|_| |_|\__,_|___/ \____/ \__\__,_|_|   \__| \_| |_/\___|_|  \___| #
#                                                                                                         #
# ------------------------------------------------------------------------------------------------------- #

# WARNING, depending on the size of database the created inputdir can be quite large.
# Estimate: for ~100 entries inputdir is about 300mb.

# args = parseArguments() # Get input arguments:
# ID = args.ID            # ID of the structure we are going to work with
# n = args.deformation    # deformation coefficient - an optional input argument
# path = '../'+ID+'/'     # Prepare folder folder where all subfolders for a single ID
# os.makedirs(path)       # will be created/executed/cleaned - working directory

warnings.filterwarnings("ignore")  # slightly dirty, but simple way to disable pymatgen complaining about lack
                                   # of element labels in POSCARs from aflowlib "UserWarning: Elements in POSCAR cannot be determined."
                                   # Warning ! This disables _ALL_ warnings.

### Setting which database we work with ###

wdatadir_structure = '../Database/datadir_structure_relaxed/'
wdatalist = '../Database/datalist_updated_sieved.mag.field_sieved.mag.sites_no.duplicates.csv'
output_path = 'D:/MCES/BEXT'

# wdatalist = '../Database/TESTS/TestDB/datalist_TestDB.csv'
# wdatadir_structure = '../Database/TESTS/TestDB/datadir_TestDB/'
# output_path = '../Database/TESTS/TestDB/datadir_TestDB/'

# wdatalist = '../Database/TESTS/TestDB_hex/datalist_TestDB_hex.csv'
# wdatadir_structure = '../Database/TESTS/TestDB_hex/datadir_TestDB_hex/'
# output_path = '../Database/TESTS/TestDB_hex/datadir_TestDB_hex/'


df = pd.read_csv(wdatalist, index_col=0, sep=',')

with tqdm.tqdm(total=len(df.index)) as pbar:  # A wrapper that creates nice progress bar
    pbar.set_description("Processing datalist")
    for item in df.index.tolist():
        pbar.update(1)
        ID = str(item)
        def_factor = 0  # deformation factor 0.05 == 5%, deformation is don both for expansion and contraction.

        path = output_path + 'inputdir/'+str(ID)+'/'     # Prepare folder folder where all subfolders for a single ID
        if os.path.exists(path):    # !WARNING! Overwrites existing path!
            shutil.rmtree(path)     #
        os.makedirs(path)           #

        ### Take all information from poscar as list of strings:
        structure_file = wdatadir_structure + str(ID)
        poscar_content = POSCAR_reader.read(structure_file)
        lat_type = (df.loc[item, 'Bravais_lattice'])
        elem_list = (df.loc[item, 'species']).replace(';', ',').replace("'", "").strip('][').split(', ')

        if def_factor == 0:            # Check if we want deformation to happen
            undeformed_lattice()       # Create a folder for a job without any deformation
        else:
            undeformed_lattice()
            deformed_lattice(lat_type, def_factor)  # Create the proper number of folders for all possible deformations according to Bravais lattice type
