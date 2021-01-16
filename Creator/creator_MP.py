import os
import shutil
import pandas as pd
import tqdm
from pymatgen.core import Lattice
from pymatgen.io.vasp import Poscar
import warnings
from Creator import INCAR_maker
from Creator import POSCAR_maker
from Creator import POTCAR_maker
from Tools import POSCAR_reader

# Placeholders for global variables
elem_list = ''
path = ''
poscar_content = ''
magmoms = ''


def single_run(calc_type, updated_lattice):
    """Creates files for a single VASP run"""

    path_to_calc = path + calc_type + "/"
    os.mkdir(path_to_calc)
    POSCAR_maker.writer_MP(path_to_calc, poscar_content, updated_lattice)
    # POSCAR_maker.structure_corrector_MP(path_to_calc)
    INCAR_maker.writer(path_to_calc, poscar_content, elem_list, calc_type, magmoms)
    POTCAR_maker.writer(path_to_calc, elem_list)


def undeformed_lattice(a, b, c, alpha, beta, gamma):
    """Creates all necessary files to run a VASP job for initial undeformed structure"""

    calculation_type = "undeformed"
    new_lattice = Lattice.from_parameters(a, b, c, alpha, beta, gamma)
    single_run(calculation_type, new_lattice)


def applied_field(a, b, c, alpha, beta, gamma):
    """Creates all necessary files to run a VASP job for structure
    under applied field (value of the field is set in INCAR_maker)"""

    calculation_type = "Applied_Field"
    new_lattice = Lattice.from_parameters(a, b, c, alpha, beta, gamma)
    single_run(calculation_type, new_lattice)


def volumetric_deformation(n, a, b, c, alpha, beta, gamma):
    deformation_types = ['V_dec', 'V_inc']
    for deformation_type in deformation_types:
        if deformation_type == "V_dec":
            new_lattice = Lattice.from_parameters(a, b, c, alpha, beta, gamma)
            volume_before = new_lattice.volume
            new_lattice = new_lattice.scale((1 - n) * volume_before)
            single_run(deformation_type, new_lattice)

        if deformation_type == "V_inc":
            new_lattice = Lattice.from_parameters(a, b, c, alpha, beta, gamma)
            volume_before = new_lattice.volume
            new_lattice = new_lattice.scale((1 + n) * volume_before)
            single_run(deformation_type, new_lattice)


def uniaxial_defromations(lattice_type, n, a, b, c, alpha, beta, gamma):
    """Takes lattice type as argument and depending on symmetry creates all necessary files to run VASP jobs for all
    possible independent deformations.
    Depending on symmetry create a 3x3 deformation matrix which determines weather a, b, c are deformed. This matrix is
    later passed to POSCAR_maker where it is used for adjusting the lattice matrix.
    Deformations are done both for increase and decrease of each lattice parameter.
    """
    # deformation factor 0.05 == 5%, deformation is done both for expansion and contraction.

    if lattice_type == "cubic":
        deformation_types = ['a_dec', 'a_inc']
        for deformation_type in deformation_types:
            if deformation_type == "a_dec":
                a = a * (1 - n)
                b = b * (1 - n)
                c = c * (1 - n)
                new_lattice = Lattice.from_parameters(a, b, c, alpha, beta, gamma)
                single_run(deformation_type, new_lattice)
            if deformation_type == "a_inc":
                a = a * (1 + n)
                b = b * (1 + n)
                c = c * (1 + n)
                new_lattice = Lattice.from_parameters(a, b, c, alpha, beta, gamma)
                single_run(deformation_type, new_lattice)

    elif lattice_type == "hexagonal":
        deformation_types = ['a_dec', 'c_dec', 'a_inc', 'c_inc']
        for deformation_type in deformation_types:
            if deformation_type == "a_dec":
                a = a * (1 - n)
                b = b * (1 - n)
                new_lattice = Lattice.from_parameters(a, b, c, alpha, beta, gamma)
                single_run(deformation_type, new_lattice)
            if deformation_type == "c_dec":
                c = c * (1 - n)
                new_lattice = Lattice.from_parameters(a, b, c, alpha, beta, gamma)
                single_run(deformation_type, new_lattice)
            if deformation_type == "a_inc":
                a = a * (1 + n)
                b = b * (1 + n)
                new_lattice = Lattice.from_parameters(a, b, c, alpha, beta, gamma)
                single_run(deformation_type, new_lattice)
            if deformation_type == "c_inc":
                c = c * (1 + n)
                new_lattice = Lattice.from_parameters(a, b, c, alpha, beta, gamma)
                single_run(deformation_type, new_lattice)

    elif lattice_type == "rhombohedral" or lattice_type == "trigonal":
        deformation_types = ['a_dec', 'a_inc']
        for deformation_type in deformation_types:
            if deformation_type == "a_dec":
                a = a * (1 - n)
                b = b * (1 - n)
                c = c * (1 - n)
                new_lattice = Lattice.from_parameters(a, b, c, alpha, beta, gamma)
                single_run(deformation_type, new_lattice)
            if deformation_type == "a_inc":
                a = a * (1 + n)
                b = b * (1 + n)
                c = c * (1 + n)
                new_lattice = Lattice.from_parameters(a, b, c, alpha, beta, gamma)
                single_run(deformation_type, new_lattice)

    elif lattice_type == "tetragonal":
        deformation_types = ['a_dec', 'c_dec', 'a_inc', 'c_inc']
        for deformation_type in deformation_types:
            if deformation_type == "a_dec":
                a = a * (1 - n)
                b = b * (1 - n)
                new_lattice = Lattice.from_parameters(a, b, c, alpha, beta, gamma)
                single_run(deformation_type, new_lattice)
            if deformation_type == "c_dec":
                c = c * (1 - n)
                new_lattice = Lattice.from_parameters(a, b, c, alpha, beta, gamma)
                single_run(deformation_type, new_lattice)
            if deformation_type == "a_inc":
                a = a * (1 + n)
                b = b * (1 + n)
                new_lattice = Lattice.from_parameters(a, b, c, alpha, beta, gamma)
                single_run(deformation_type, new_lattice)
            if deformation_type == "c_inc":
                c = c * (1 + n)
                new_lattice = Lattice.from_parameters(a, b, c, alpha, beta, gamma)
                single_run(deformation_type, new_lattice)

    elif lattice_type == "orthorhombic":
        deformation_types = ['a_dec', 'b_dec', 'c_dec', 'a_inc', 'b_inc', 'c_inc']
        for deformation_type in deformation_types:
            if deformation_type == "a_dec":
                a = a * (1 - n)
                new_lattice = Lattice.from_parameters(a, b, c, alpha, beta, gamma)
                single_run(deformation_type, new_lattice)
            if deformation_type == "b_dec":
                b = b * (1 - n)
                new_lattice = Lattice.from_parameters(a, b, c, alpha, beta, gamma)
                single_run(deformation_type, new_lattice)
            if deformation_type == "c_dec":
                c = c * (1 - n)
                new_lattice = Lattice.from_parameters(a, b, c, alpha, beta, gamma)
                single_run(deformation_type, new_lattice)
            if deformation_type == "a_inc":
                a = a * (1 + n)
                new_lattice = Lattice.from_parameters(a, b, c, alpha, beta, gamma)
                single_run(deformation_type, new_lattice)
            if deformation_type == "b_inc":
                b = b * (1 + n)
                new_lattice = Lattice.from_parameters(a, b, c, alpha, beta, gamma)
                single_run(deformation_type, new_lattice)
            if deformation_type == "c_inc":
                c = c * (1 + n)
                new_lattice = Lattice.from_parameters(a, b, c, alpha, beta, gamma)
                single_run(deformation_type, new_lattice)

    elif lattice_type == "monoclinic":
        deformation_types = ['a_dec', 'b_dec', 'c_dec', 'a_inc', 'b_inc', 'c_inc']
        for deformation_type in deformation_types:
            if deformation_type == "a_dec":
                a = a * (1 - n)
                new_lattice = Lattice.from_parameters(a, b, c, alpha, beta, gamma)
                single_run(deformation_type, new_lattice)
            if deformation_type == "b_dec":
                b = b * (1 - n)
                new_lattice = Lattice.from_parameters(a, b, c, alpha, beta, gamma)
                single_run(deformation_type, new_lattice)
            if deformation_type == "c_dec":
                c = c * (1 - n)
                new_lattice = Lattice.from_parameters(a, b, c, alpha, beta, gamma)
                single_run(deformation_type, new_lattice)
            if deformation_type == "a_inc":
                a = a * (1 + n)
                new_lattice = Lattice.from_parameters(a, b, c, alpha, beta, gamma)
                single_run(deformation_type, new_lattice)
            if deformation_type == "b_inc":
                b = b * (1 + n)
                new_lattice = Lattice.from_parameters(a, b, c, alpha, beta, gamma)
                single_run(deformation_type, new_lattice)
            if deformation_type == "c_inc":
                c = c * (1 + n)
                new_lattice = Lattice.from_parameters(a, b, c, alpha, beta, gamma)
                single_run(deformation_type, new_lattice)

    elif lattice_type == "triclinic":
        deformation_types = ['a_dec', 'b_dec', 'c_dec', 'a_inc', 'b_inc', 'c_inc']
        for deformation_type in deformation_types:
            if deformation_type == "a_dec":
                a = a * (1 - n)
                new_lattice = Lattice.from_parameters(a, b, c, alpha, beta, gamma)
                single_run(deformation_type, new_lattice)
            if deformation_type == "b_dec":
                b = b * (1 - n)
                new_lattice = Lattice.from_parameters(a, b, c, alpha, beta, gamma)
                single_run(deformation_type, new_lattice)
            if deformation_type == "c_dec":
                c = c * (1 - n)
                new_lattice = Lattice.from_parameters(a, b, c, alpha, beta, gamma)
                single_run(deformation_type, new_lattice)
            if deformation_type == "a_inc":
                a = a * (1 + n)
                new_lattice = Lattice.from_parameters(a, b, c, alpha, beta, gamma)
                single_run(deformation_type, new_lattice)
            if deformation_type == "b_inc":
                b = b * (1 + n)
                new_lattice = Lattice.from_parameters(a, b, c, alpha, beta, gamma)
                single_run(deformation_type, new_lattice)
            if deformation_type == "c_inc":
                c = c * (1 + n)
                new_lattice = Lattice.from_parameters(a, b, c, alpha, beta, gamma)
                single_run(deformation_type, new_lattice)


def creator(datalist, datadir, def_factor=0.05, undeformed=True, hydrostatic=False, uniaxial=False, field=False):
    df = pd.read_csv(datalist, index_col=0, sep=',')
    with tqdm.tqdm(total=len(df.index)) as pbar:  # A wrapper that creates nice progress bar
        pbar.set_description("Processing datalist")
        for item in df.index.tolist():
            pbar.update(1)
            ID = str(item)
            global path
            path = output_path + '/' + str(ID) + '/'  # Prepare folder single ID (all subfolders will go there)
            if os.path.exists(path):  # !WARNING! Overwrites existing path!
                shutil.rmtree(path)  #
            os.makedirs(path)  #

            # Take all information from poscar as list of strings:

            structure_file = datadir + str(ID) + '.cif'
            global poscar_content
            poscar_content = POSCAR_reader.read(structure_file)
            poscar_string = ''.join(poscar_content)  # merging poscar content in a single string so pymatgen can read it
            poscar = Poscar.from_string(poscar_string)  # using pymatgen to acquire our structure from poscar content
            structure = poscar.structure


            a_par = structure.lattice.a
            b_par = structure.lattice.b
            c_par = structure.lattice.c
            alpha_ang = structure.lattice.alpha
            beta_ang = structure.lattice.beta
            gamma_ang = structure.lattice.gamma

            lat_type = (df.loc[item, 'lattice_system'])

            # getting list of elemets from structure
            global elem_list
            used = set()
            mylist = structure.species
            unique = [x for x in mylist if x not in used and (used.add(x) or True)]

            elem_list = str(unique).replace('Element ', '').strip('][').split(', ')


            # global magmoms
            # magmoms = (df.loc[item, 'magmom']).replace(';', '').strip('][')

            if undeformed is True:  # check if we want undeformed structure
                undeformed_lattice(a_par, b_par, c_par, alpha_ang, beta_ang, gamma_ang)
            if uniaxial is True & (def_factor != 0):  # Check if we want deformation to happen
                uniaxial_defromations(lat_type, def_factor, a_par, b_par, c_par, alpha_ang, beta_ang, gamma_ang)
            if (hydrostatic is True) & (def_factor != 0):  # check if we want uniaxial deformation
                volumetric_deformation(def_factor, a_par, b_par, c_par, alpha_ang, beta_ang, gamma_ang)
            if field is True:  # check if we want applied field on undeformed structure
                applied_field(a_par, b_par, c_par, alpha_ang, beta_ang, gamma_ang)


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
                                   # of element labels in POSCARs from aflowlib
                                   # "UserWarning: Elements in POSCAR cannot be determined."
                                   # Warning ! This disables _ALL_ warnings.

### Setting which database we work with ###

wdatadir = 'D:/MCES/MP/datadir/'
# wdatalist = 'D:/MCES/MP/step4_failed_sieved_input_for_step5.csv'

wdatalist = 'D:/MCES/MP/step0.csv'

# wdatalist = '../Database/MP/dt.csv'
# wdatalist = 'D:/MCES/MP/datalist_lattfix_updated_sieved.mag.field_sieved.mag.sites_no.duplicates.csv'


output_path = 'D:/MCES/MP/inputdir_V'


# wdatadir = 'D:/MCES/MP/test4/MP_structures/'
# wdatalist = 'D:/MCES/MP/test4/datalist_MP_lattfix.csv'
# output_path = 'D:/MCES/MP/test4/'

# wdatalist = '../Database/TESTS/TestDB/datalist_TestDB.csv'
# wdatadir_structure = '../Database/TESTS/TestDB/datadir_TestDB/'
# output_path = '../Database/TESTS/TestDB/datadir_TestDB/'

# wdatalist = '../Database/TESTS/TestDB_hex/datalist_TestDB_hex.csv'
# wdatadir_structure = '../Database/TESTS/TestDB_hex/datadir_TestDB_hex/'
# output_path = '../Database/TESTS/TestDB_hex/datadir_TestDB_hex/'

creator(wdatalist, wdatadir, def_factor=0.05, hydrostatic=True, undeformed=False, field=False)
# creator('MnAs.csv', 'datadir/')