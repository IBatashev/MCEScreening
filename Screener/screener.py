import pandas as pd
import phonopy
from ase import io

from tools import POSCAR_reader


def calculate_mag_field(moment, volume):
    """Takes cell moment in [mB] and volume in [A^3] and returns value for internal magnetic field in [T]"""

    pi = 3.141592653
    mu0 = 4*pi*10**(-7)                        # Vacuum permeability in H/m
    mB = 9.2741*10**(-24)                      # Bohr magneton value in J/T
    field = (mu0*moment*mB)/(volume*10**(-30)) # Formula for internal magnetic field in Tesla
    return field

def check_universal_scaling_factor(ID):
    """Checks if universal scaling factor constant was modified in POSCAR (!= 1.0) for the structure
     and returns either True + value of scaling factor or False"""

    poscar_file = '../Database/poscars_for_tests/'+str(ID) # for tests CHANGE LATER

    poscar_file = '../Database/sample_datadir/'+str(ID)
    poscar_content = POSCAR_reader.read(poscar_file)
    scaling_factor = float(poscar_content[1])              # read universal scaling factor from poscar (second line in all POSCARs)
    if scaling_factor != 1.0:                              # Check if equals 1.0  - perhaps I can add some tolerance here...
        return True, scaling_factor                        # Return True if modified and value of universal scaling factor
    else:
        return False, ''

def mag_sites_calculator(ID):
    """Determines how many unique magnetic sites are present in the structure.
    Takes ID as input, looks up poscar file in the datadir and returns number of sites as integer"""

    # poscar_file = '../Database/datadir/'+str(ID)
    poscar_file = '../Database/poscars_for_tests/'+str(ID) # for tests CHANGE LATER

    poscar_file = '../Database/sample_datadir/'+str(ID)

    ### List of Magnetic Atoms
    # need to add all of them later
    magnetic = ['Cr', 'Mn', 'Fe', 'Co', 'Ni', 'Gd']

    ### Open POSCAR with phonopy and get Wyckoff letters and chemical symbols from it
    # maybe should try with ase????????
    ph = phonopy.load(unitcell_filename=poscar_file)
    sym_list = ph.get_symmetry().get_Wyckoff_letters()
    at_type_len = len(ph.unitcell.get_positions())
    at_type_list = []

    ### Open POSCAR and get element symbols
    poscar_content = POSCAR_reader.read(poscar_file)

    for i in range(7, 7 + at_type_len):             # POSCAR files from aflowlib have 7 lines of text before atomic coordinates! so we start at line 8
        l = str.split(poscar_content[i])            # Aflow POSCAR has symbols for elements listed after their coordinate
        at_type_list.append(str(l[3]))              # create a list of all elements in the structure

    ### Check how many unique MAGNETIC sites are present
    at_list = []
    for num, val in enumerate(sym_list):            # build list of magnetic positions
        if at_type_list[num] in magnetic:
           at_list.append(val)
    num_unique = len(set(at_list))                  # check how many are unique
    return num_unique                               # return resulting number of sites


    ### Another attempt to do same as above need to test both further

    # print(ph.get_symmetry().get_dataset())
    # for i, vall in enumerate(at_type_list):
    #     sym_point = ph.get_symmetry().get_site_point_group()
    #     # sym_list = sym_list.append(ph.get_symmetry().get_pointgroup(sym_point))
    #     print(sym_point)
    # print(sym_list)
    # print(at_type_list)

def mag_sites_difference(ID):
    # Works after calculation, reads OUTCAR and checks values of moments for different magnetic sites and returns biggest
    # difference between moments of sublattices or just all significantly different moments...
    return

def duplicates(ID):
    # report how many different lattice types are there for the chem formula => may indicate instability
    # for same structure only leave one, but how we decide...
    # Also sometimes reports for structure types are same just lower level symmetry... how do we deal with it? maybe
    # compare atomic positions list...
    return

def screener_before(datalist):
    """First main function that works with screening database. Used to apply criteria thet don't require calculations.
    Loops over all database entries and applies initial screening parameters/checks.
    Creates a shorter database that we will use for submitting calculations"""

    df = pd.read_csv(datalist, index_col=0, sep=',')
    for item in df.index.tolist():
        ### Write number of sites into datalist:
        df.loc[item, 'mag_sites'] = mag_sites_calculator(item)
        ### Write internal magnetic field into datalist:
        df.loc[item, 'mag_field'] = calculate_mag_field(df.loc[item, 'moment_cell'], df.loc[item, 'volume_cell'])
        ### Check if universal scaling factor is 1.0, othervise results for magnetic field calculaded using aflow data are unreliable
        is_scalled, scaling_factor = check_universal_scaling_factor(item)
        if is_scalled == True:
            df.loc[item, 'comment1'] = 'scaling factor = ' + str(scaling_factor)
        ### Work with duplicates...

        ### Write to a separate shorter file
        df.to_csv(datalist.replace(".csv",'_edited.csv'))
        ### copies corresponding POSCAR files to a new smaller datadir tham we will move to cluster

        # for num, val in enumerate(sym_list):
        # if at_type_list[num] in Rb:
        #     df = df.drop([item], axis=0)

def screener_after(datalist):
    # Second main function used for applying criteria obtained after calculations ad performing corresponding checks
    df = pd.read_csv(datalist, index_col=0, sep=',')
    for item in df.index.tolist():
        mag_sites_difference(item)
        #...

screener_before('../Database/sample_datalist.csv')
# print(mag_sites_calculator(1770))
# print(check_universal_scaling_factor(24))
