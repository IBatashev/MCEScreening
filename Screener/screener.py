import pandas as pd
import phonopy
from ase import io
import shutil
import os
import numpy as np
### LOCAL IMPORTS ###
from tools import POSCAR_reader

### Setting which database we work with ###
wdatadir_structure = '../Database/datadir_structure/'
wdatadir_aflow = '../Database/datadir_aflow/'
wdatadir_incars = '../Database/datadir_incars/'
wdatalist = '../Database/datalist.csv'

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

    structure_file = wdatadir_structure+str(ID)
    poscar_content = POSCAR_reader.read(structure_file)
    scaling_factor = float(poscar_content[1])              # read universal scaling factor from poscar (second line in all POSCARs)
    if scaling_factor != 1.0:                              # Check if equals 1.0  - perhaps I can add some tolerance here...
        return True, scaling_factor                        # Return True if modified and value of universal scaling factor
    else:
        return False, ''

def mag_sites_calculator(ID):
    """Determines how many unique magnetic sites are present in the structure.
    Takes ID as input, looks up structure file in the datadir and returns number of unique sites as integer"""

    ### List of Magnetic Atoms
    # need to add all of them later
    # magnetic = ['Cu','Cr', 'Mn', 'Fe', 'Co', 'Ni', 'Gd']
    magnetic = [ 'Sc', 'Ti', 'V', 'Cr', 'Mn', 'Fe', 'Co', 'Ni', 'Cu', 'Zn', # Sc, Y, and everythong else is considered not magnetic
                 'Y', 'Zr', 'Nb', 'Mo', 'Tc', 'Ru', 'Rh', 'Pd', 'Ag', 'Cd', # we sieve out Cd at earliear step, so I could have ommitted it from this list but this approach is more "universal"
                 'La', 'Ce', 'Pr', 'Nd', 'Pm', 'Sm', 'Eu', 'Gd', 'Tb', 'Dy', 'Ho', 'Er', 'Tm', 'Yb', 'Lu']
    bad_notation = False                                                  # sometimes aflow files have alphabet instead of element symbols in structure POSCAR I call it bad notation
    rename_list = ['A', 'B', 'C', 'D', 'E', 'F', 'G', 'H', 'I', 'J', 'K'] # A list for fixing bad notation
                                                                          # I don't expect there to be a composition with more than 11 components so here I stop at K
    elemlist = []
    site_counter = 0

    df = pd.read_csv(wdatalist, index_col=0, sep=',')

    #### NEW APPROACH ###
    structure_file = wdatadir_structure + str(ID)         # still need to decide weather I should
                                                          # use relaxed or non-relaxed structure for this - non relaxed may give overestimaded positions, if symmetry was lowered hz
    with open(structure_file, 'r') as f:    # need to get all atoms from Wyckoff subsection of structure file
        for line in f:
            if 'Representative' in line:    # subsection we are interested in goes after line with word "Representative"
                for line in f:              # now you are at the lines you want
                    if 'WYCCAR' in line:    # Ends before line with "WYCCAR"
                        break
                    else:
                        elem = line.split()[3]  # element symbol is after Wyckoff x y z coordinates so 4th item in list
                        if elem == 'A':         # check if notation is shitty - in this case first element is represented with letter A (fortunately there is no element in periodic table labeled with A)
                            bad_notation = True
                        if bad_notation == True:# Apply correction for bad notation replacing all alphabet letters with proper chemical symbols taken from list of elements in composition
                            n = rename_list.index(elem)
                            temp_line = (df.loc[ID, 'species']).strip('[').strip(']').strip(' ')
                            temp_line2 = temp_line.replace("'", "")
                            species_list = np.asarray(temp_line2.split('; '))
                            elem = species_list[n]
                        elemlist = np.append(elemlist, elem)
    for k in elemlist:
        if k in magnetic:
            site_counter = site_counter + 1
    return site_counter

    #####################

    # poscar_file = '../Database/datadir/'+str(ID)
    # poscar_file = '../Database/poscars_for_tests/'+str(ID) # for tests CHANGE LATER
    #
    # poscar_file = '../Database/sample_datadir/'+str(ID)

    # ### Open POSCAR with phonopy and get Wyckoff letters and chemical symbols from it
    # # maybe should try with ase????????
    # ph = phonopy.load(unitcell_filename=poscar_file)
    # sym_list = ph.get_symmetry().get_Wyckoff_letters()
    # at_type_len = len(ph.unitcell.get_positions())
    # at_type_list = []
    #
    # ### Open POSCAR and get element symbols
    # poscar_content = POSCAR_reader.read(poscar_file)
    #
    # for i in range(7, 7 + at_type_len):             # POSCAR files from aflowlib have 7 lines of text before atomic coordinates! so we start at line 8
    #     l = str.split(poscar_content[i])            # Aflow POSCAR has symbols for elements listed after their coordinate
    #     at_type_list.append(str(l[3]))              # create a list of all elements in the structure
    #
    # ### Check how many unique MAGNETIC sites are present
    # at_list = []
    # for num, val in enumerate(sym_list):            # build list of magnetic positions
    #     if at_type_list[num] in magnetic:
    #        at_list.append(val)
    # num_unique = len(set(at_list))                  # check how many are unique
    # return num_unique                               # return resulting number of sites

    #########################################

    ### Another attempt to do same as above need to test both further ###

    # print(ph.get_symmetry().get_dataset())
    # for i, vall in enumerate(at_type_list):
    #     sym_point = ph.get_symmetry().get_site_point_group()
    #     # sym_list = sym_list.append(ph.get_symmetry().get_pointgroup(sym_point))
    #     print(sym_point)
    # print(sym_list)
    # print(at_type_list)

    #########################################

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
    # I feel that performance may not be optimal - making two loops does not seem reasonable
    # but I have to test if working with two df simultaneously is faster...
    """First main function that works with screening database. Used to apply criteria thet don't require calculations.
    Loops over all database entries and applies initial screening parameters/checks.
    Creates a shorter database that we will use for submitting calculations"""

    # Criteria limits
    min_site_number = 1
    min_mag_field = 0.45

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

        print(item)
     ### Wrire an updated datalist back to file
    df.to_csv(datalist)


    # ### Write the shorter sieved database to a separate file and copy relevant POSCAR to new datadir
    # if os.path.exists('../Database/datadir_sieved/'):  # Prepare new datadir folder
    #     shutil.rmtree('../Database/datadir_sieved/')  # (cleans old one if it already existed)
    # os.makedirs('../Database/datadir_sieved/')

    ### loops through database and drops everything that does not fit the criteria
    dropper(datalist, 'mag_filed', min_mag_field)
    # dropper(datalist, 'mag_sites', min_site_number)

def dropper(datalist, sieve_type, sieve_size):
    """Removes entries from datalist according to selected sieving ctriteria"""

    df = pd.read_csv(datalist, index_col=0, sep=',')
    for item in df.index.tolist():
        if df.loc[item, sieve_type] <= sieve_size:
            df = df.drop([item], axis=0)
        # else: # copies corresponding  files to a new smaller datadir that we will move to cluster
        #     shutil.copy('../Database/datadir/'+str(item), '../Database/datadir_sieved/'+str(item))
    df.to_csv(datalist.replace(".csv", '_sieved' + sieve_type + '.csv'))

def screener_after(datalist):
    # Second main function used for applying criteria obtained after calculations ad performing corresponding checks
    df = pd.read_csv(datalist, index_col=0, sep=',')
    for item in df.index.tolist():
        mag_sites_difference(item)
        #...mag field again - to correct for erroneously scalled compounds


### TEMP WORKING AREA ###
# df = pd.read_csv('../Database/datalist.csv', index_col=0, sep=',')
# # f = open('site_count_list.txt', 'w+')
# table = []
# for item in df.index.tolist():
#     table = np.append(table, mag_sites_calculator(item))
#     print(item)
# np.savetxt('site_count_list_relaxed.txt', table, delimiter='\n')

# print(mag_sites_calculator(0)) #1938
# print(check_universal_scaling_factor(24))

# screener_before(wdatalist)

print(POSCAR_reader.read(wdatadir_structure+'0'))