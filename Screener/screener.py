import pandas as pd
import phonopy
from ase import io

def screener(datalist):
    df = pd.read_csv(datalist, index_col=0, sep=',')
    for item in df.index.tolist():
        ### Write number of sites into datalist:
        df.loc[item, 'mag_sites'] = mag_sites_calculator(item['ID'])
        ### Write number of sites into datalist:
        df.loc[item,'mag_field'] = calculate_mag_field(df.loc[item,'moment_cell'], df.loc[item,'volume_cell'])
        ### Work with duplicates...

        
    # for num, val in enumerate(sym_list):
    # if at_type_list[num] in Rb:
    #     df = df.drop([item], axis=0)

    # df.to_csv()

def calculate_mag_field(moment, volume):
    """Takes cell moment in [mB] and volume in [A^3] and returns value for internal magnetic field in [T]"""

    pi = 3.141592653
    mu0 = 4*pi*10**(-7)                        # Vacuum permeability in H/m
    mB = 9.2741*10**(-24)                      # Bohr magneton value in J/T
    field = (mu0*moment*mB)/(volume*10**(-30)) # Formula for internal magnetic field in Tesla
    return field


def mag_sites_calculator(ID):
    """Determines how many unique magnetic sites are present in the structure.
    Takes ID as input, looks up poscar file in the datadir and returns number of sites as integer"""

    poscar_file = '../Database/datadir/'+str(ID)
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
    # perhaps make this a separrate function for both runner and screener...
    poscar = open(poscar_file, "r")
    poscar_content = []
    for line in poscar:
        poscar_content.append(line)
    poscar.close()

    for i in range(7, 7 + at_type_len):             # POSCAR files from aflowlib have 7 lines of text before atomic coordinates! so we start at line 8
        l = str.split(poscar_content[i], '  ')      # Aflow POSCAR has symbols for elements listed after their coordinate
        at_type_list.append(str(l[4]))              # create a list of all elements in the structure

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