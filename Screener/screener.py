import pandas as pd
import shutil
import os
import numpy as np
import tqdm
import tarfile
### LOCAL IMPORTS ###
import POSCAR_reader


### Setting which database we work with ###
wdatadir_structure = '../Database/datadir_structure_relaxed/'
wdatadir_aflow = '../Database/datadir_aflow/'
wdatadir_incars = '../Database/datadir_incars/'
# wdatalist = '../Database/datalist.csv'
# wdatalist = '../Database/datalist_updated_sieved.mag.field_sieved.mag.sites.csv'
vasp_results_dir = '../Database/TestDB/second_run/1.2/outdir'
wdatalist = '../Database/TestDB/datalist_TestDB.csv'


def calculate_mag_field(moment, volume):
    """Takes cell moment in [mB] and volume in [A^3] and returns value for internal magnetic field in [T]"""

    pi = 3.141592653
    mu0 = 4*pi*10**(-7)                         # Vacuum permeability in H/m
    mB = 9.2741*10**(-24)                       # Bohr magneton value in J/T
    field = (mu0*moment*mB)/(volume*10**(-30))  # Formula for internal magnetic field in Tesla
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
    magnetic = [ 'Ti', 'V', 'Cr', 'Mn', 'Fe', 'Co', 'Ni', 'Cu',  # Sc, Y, and everythong else is considered not magnetic
                 'Nb', 'Mo', 'Ru', 'Rh', 'Pd', # we sieve out Cd at earliear step, so I could have ommitted it from this list but this approach is more "universal"
                 'Ce', 'Pr', 'Nd', 'Sm', 'Eu', 'Gd', 'Tb', 'Dy', 'Ho', 'Er', 'Tm', 'Yb']
    bad_notation = False                                                    # sometimes aflow files have alphabet instead of element symbols in structure POSCAR I call it bad notation
    rename_list = ['A', 'B', 'C', 'D', 'E', 'F', 'G', 'H', 'I', 'J', 'K']   # A list for fixing bad notation
                                                                            # I don't expect there to be a composition with more than 11 distinctive components so here I stop at K

    eltable = np.empty([0, 2])
    df = pd.read_csv(wdatalist, index_col=0, sep=',')

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
                        wyckoff = line.split()[4]+line.split()[5]
                        if elem == 'A':         # check if notation is shitty - in this case first element is represented with letter A (fortunately there is no element in periodic table labeled with A)
                            bad_notation = True
                        if bad_notation == True:# Apply correction for bad notation replacing all alphabet letters with proper chemical symbols taken from list of elements in composition
                            n = rename_list.index(elem)
                            temp_line = (df.loc[ID, 'species']).strip('[').strip(']').strip(' ')
                            temp_line2 = temp_line.replace("'", "")
                            species_list = np.asarray(temp_line2.split('; '))
                            elem = species_list[n]
                        eltable = np.append(eltable, [[elem, wyckoff]], axis=0)
    # This section with masks probably could be done better and shorter but with such smal arrays it shouldn't matter much
    mask = np.in1d(eltable[:, 0], magnetic)  # Check what entries in Wyckoff list are in list of magnetic atoms
    eltable2 = eltable[mask]                 # new array of only magnetic ones using previous mask
    unique_keys, mask2 = np.unique(eltable2[:, 1], return_index=True)  # check what entries of the new array have unique wyckoff sites - creates a list of indicies corresponding to said entries
    site_counter = np.size(eltable2[mask2], axis=0)                    # get number of unique sites counting number of obtained indicies
    return site_counter


def duplicates(datalist):
    """Remove duplicates from datalist. Entries are considered duplicates if they have the same composition AND same
     spacegroup. Duplicates with different spacegroups are kept - a new field is added to datalist counting how many
     structural variations are present for same composition (May serve as indication of instability)."""

    # Sometimes reports for structure types are same just lower level symmetry... how do we deal with it? maybe
    # compare atomic positions list...

    # I will not remove items with doubled etc chemical structure e.g. both Fe2O and Fe4O2 will remain as in some cases
    # in aflow new wyckoff sites were attributed to these extra atoms and I am not confident in judging how significant
    # they are. Perhaps a more thorough check? For now it doesn't matter but for bigger database it may become important
    # to reduce sample size as much as possible

    df = pd.read_csv(datalist, index_col=0, sep=',')
    df = df.drop_duplicates(subset=['compound', 'spacegroup'])  # remove all entries with exact same spacegroup and composition - their difference is most likely a sligh change of lattice parameters
    df['polymorphs'] = (df.groupby('compound').compound.transform(lambda x: x.duplicated(keep=False).sum()))  # check how many entries still have the same composition - are polymorphs with diffeent structures - and count them
    df.to_csv(datalist.replace('.csv', '_no.duplicates.csv'))


def timing(ID):
    """Reads OUTCARs and returns time taken to calculate entry (total: undeformed + all deformations)"""

    time = 0
    deformations = os.listdir(vasp_results_dir + '/' + str(ID))
    for deformation in deformations:
        with open(vasp_results_dir + '/' + str(ID) + '/' + deformation + "/OUTCAR") as search:
            for line in search:
                if 'User time' in line:
                    time = time + (float(line.split()[-1]))
    return round(time)

def memmory_used(ID):
    """Reads OUTCARs and returns memmory used to calculate entry (total: undeformed + all deformations)"""

    memmory = 0
    deformations = os.listdir(vasp_results_dir + '/' + str(ID))
    for deformation in deformations:
        with open(vasp_results_dir + '/' + str(ID) + '/' + deformation + "/OUTCAR") as search:
            for line in search:
                if 'Maximum memory used' in line:
                    memmory = memmory + (float(line.split()[-1]))
    memmory = memmory * 9.5367432e-7  # converting kb to Gb
    return round(memmory, 2)


def geometry_after(ID, deformation):
    """Reads OUTCARs and returns lattice parameters """

    if deformation in os.listdir(vasp_results_dir + '/' + str(ID)):
        with open(vasp_results_dir + '/' + str(ID) + '/' + deformation + "/OUTCAR") as search:
            a = 0
            b = 0
            c = 0
            for line in search:
                if 'ALAT' in line:
                    a = round((float(line.split()[-1])), 2)
                if 'C/A-ratio' in line:
                    c = round((float(line.split()[-1]))*a, 2)
                if 'B/A-ratio' in line:
                    b = round((float(line.split()[-1]))*a, 2)
        return a, b, c
    else:
        return 0, 0, 0


def moment_volume_after(ID, deformation):
    """Read OUTCAR file for chosen deformation type and return value of magnetic field based on moment and volume.
    If chosen deformation is n/a for this ID a '-' string will be returned"""

    if deformation in os.listdir(vasp_results_dir + '/' + str(ID)):
        with open(vasp_results_dir + '/' + str(ID) + '/' + deformation + "/OUTCAR") as search:
            moment_lines = []
            for line in search:
                if 'volume of cell' in line:
                    volume = round((float(line.split()[-1])), 2)
                if 'tot' in line:
                    moment_lines.append(line)
        moment = round(float(moment_lines[-2].split()[-1]), 2)
        field_after = calculate_mag_field(moment, volume)
        return moment, volume, field_after
    else:
        return 0, 0, 0


def sym_after(ID, deformation):
    if deformation in os.listdir(vasp_results_dir + '/' + str(ID)):
        with open(vasp_results_dir + '/' + str(ID) + '/' + deformation + "/OUTCAR") as search:
            for line in search:
                if 'Routine SETGRP: Setting up the symmetry group for a' in line:
                    symmetry = (next(search)).strip("\n")
        return symmetry
    else:
        return 0


def archive_reader(ID, deformation):
    if deformation in os.listdir(vasp_results_dir + '/' + str(ID)):
        archive_path = vasp_results_dir + '/' + str(ID) + '/' + deformation + '/'
        for item in os.listdir(archive_path):
            if '.tar' in item:
                tar_file = item
        tar = tarfile.open(archive_path + tar_file, "r:gz")
        for member in tar.getmembers():
            if 'OSZICAR' in str(member):
                f = tar.extractfile(member)
                if f is not None:
                    content = str(f.read()).split('\\n')
                    algorithm_step = content[-3].split()
        return algorithm_step[0] + algorithm_step[1]
    else:
        return 0


def energy(ID, deformation):
    if deformation in os.listdir(vasp_results_dir + '/' + str(ID)):
        with open(vasp_results_dir + '/' + str(ID) + '/' + deformation + "/OUTCAR") as search:
            energy_lines = []
            for line in search:
                if 'free  energy   TOTEN' in line:
                    energy_lines.append(line)
        tot_energy = round(float(energy_lines[-1].split()[-2]), 4)
        return tot_energy
    else:
        return 0


def mag_sites_difference(ID):
    # Works after calculation, reads OUTCAR and checks values of moments for different magnetic sites and returns biggest
    # difference between moments of sublattices or just all significantly different moments...
    return


def screener_before(datalist):
    # I feel that performance may not be optimal - making two loops does not seem reasonable
    # but I have to test if working with two df simultaneously is faster...
    """First main function that works with screening database. Used to apply criteria that don't require calculations.
    Loops over all database entries and applies initial screening parameters/checks.
    Creates a shorter database that we will use for submitting calculations"""

    # Criteria limits
    min_site_number = 1
    min_mag_field = 0.45

    df = pd.read_csv(datalist, index_col=0, sep=',')
    with tqdm.tqdm(total=len(df.index)) as pbar:        # A wrapper that creates nice progress bar
        pbar.set_description("Processing datalist")
        for item in df.index.tolist():
            pbar.update(1)                              # Updating progress bar at each step
        ### Write number of sites into datalist:
            df.loc[item, 'mag_sites'] = mag_sites_calculator(item)
        ### Write internal magnetic field into datalist:
            df.loc[item, 'mag_field'] = calculate_mag_field(df.loc[item, 'moment_cell'], df.loc[item, 'volume_cell'])
        ### Check if universal scaling factor is 1.0, othervise results for magnetic field calculaded using aflow data are unreliable
            is_scalled, scaling_factor = check_universal_scaling_factor(item)
            if is_scalled == True:
                df.loc[item, 'comment1'] = 'scaling factor = ' + str(scaling_factor)
        ### Wrire an updated datalist back to file
        df.to_csv(datalist.replace('.csv', '_updated.csv'))

        # ### Write the shorter sieved database to a separate file and copy relevant POSCAR to new datadir
        # if os.path.exists('../Database/datadir_sieved/'):  # Prepare new datadir folder
        #     shutil.rmtree('../Database/datadir_sieved/')  # (cleans old one if it already existed)
        # os.makedirs('../Database/datadir_sieved/')

        ### drops everything that does not fit the criteria, creating new datalist file at each sieve
    sieve(datalist.replace('.csv', '_updated.csv'), 'mag_field', min_mag_field)
    sieve(datalist.replace('.csv', '_updated_sieved.mag.field.csv'), 'mag_sites', min_site_number)
        ### drops duplicates, creates new datalist file __ PERHAPS THIS SHOULD BE FIRST STET TO LOWER FURTHER WORKLOAD
    duplicates(datalist.replace('.csv', '_updated_sieved.mag.field_sieved.mag.sites.csv'))
    print("Done")


def sieve(datalist, sieve_type, sieve_size):
    """Removes entries from datalist according to selected sieving criteria"""

    print("\nSieving by " + str(sieve_type) + " with cutoff set as " + str(sieve_size))
    df = pd.read_csv(datalist, index_col=0, sep=',')
    for item in df.index.tolist():
        if df.loc[item, sieve_type] <= sieve_size:
            df = df.drop([item], axis=0)
    df.to_csv(datalist.replace(".csv", '_sieved.' + sieve_type.replace('_', '.') + '.csv'))


def screener_after(datalist):
    """Second main function used for applying criteria obtained after calculations and performing corresponding checks"""

    df = pd.read_csv(datalist, index_col=0, sep=',')
    with tqdm.tqdm(total=len(df.index)) as pbar:        # A wrapper that creates nice progress bar
        pbar.set_description("Processing datalist")
        for item in df.index.tolist():
            pbar.update(1)                              # Updating progress bar at each step

        ### Writes memmory used for calculation
            df.loc[item, 'memmory_used'] = memmory_used(item)
        ### Writes time used for calculation
            df.loc[item, 'time_to_calculate'] = timing(item)
        ### Writes magnetic fields for undeformed and deformed structures from calculation results
            #df.loc[item, 'magF_diff_a'] = (calculate_mag_field(moment_volume_after(item, 'a')[0] - moment_volume_after(item, 'undeformed')[0], moment_volume_after(item, 'a')[1] - moment_volume_after(item, 'undeformed')[1]))
            df.loc[item, 'moment_u'] = moment_volume_after(item, 'undeformed')[0]
            df.loc[item, 'moment_a'] = moment_volume_after(item, 'a')[0]
            df.loc[item, 'moment_b'] = moment_volume_after(item, 'b')[0]
            df.loc[item, 'moment_c'] = moment_volume_after(item, 'c')[0]
            df.loc[item, 'moment_V'] = moment_volume_after(item, 'V')[0]
            #
            df.loc[item, 'volume_u'] = moment_volume_after(item, 'undeformed')[1]
            df.loc[item, 'volume_a'] = moment_volume_after(item, 'a')[1]
            df.loc[item, 'volume_b'] = moment_volume_after(item, 'b')[1]
            df.loc[item, 'volume_c'] = moment_volume_after(item, 'c')[1]
            df.loc[item, 'volume_V'] = moment_volume_after(item, 'V')[1]
            #
            df.loc[item, 'a_u'] = geometry_after(item, 'undeformed')[0]
            df.loc[item, 'b_u'] = geometry_after(item, 'undeformed')[1]
            df.loc[item, 'c_u'] = geometry_after(item, 'undeformed')[2]

            df.loc[item, 'a_a'] = geometry_after(item, 'a')[0]
            df.loc[item, 'b_b'] = geometry_after(item, 'b')[1]
            df.loc[item, 'c_c'] = geometry_after(item, 'c')[2]
            #
            df.loc[item, 'magF_u'] = moment_volume_after(item, 'undeformed')[2]
            df.loc[item, 'magF_a'] = moment_volume_after(item, 'a')[2]
            df.loc[item, 'magF_b'] = moment_volume_after(item, 'b')[2]
            df.loc[item, 'magF_c'] = moment_volume_after(item, 'c')[2]
            df.loc[item, 'magF_V'] = moment_volume_after(item, 'V')[2]
        ### placeholder for futureimplementation of site difference check
            mag_sites_difference(item)
            df.loc[item, 'sym_after'] = sym_after(item, 'undeformed')
            df.loc[item, 'sym_a'] = sym_after(item, 'a')
            df.loc[item, 'sym_b'] = sym_after(item, 'b')
            df.loc[item, 'sym_c'] = sym_after(item, 'c')
            df.loc[item, 'sym_V'] = sym_after(item, 'V')
            # Get energy
            df.loc[item, 'energy'] = energy(item, 'undeformed')
            # get algorithm step at which optimization stopped
            df.loc[item, 'step'] = archive_reader(item, 'undeformed')
    df.to_csv(datalist.replace(".csv", '_out' + '.csv'))


###  WORKING AREA ###
wdatalist = '../Database/TestDB/datalist_TestDB.csv'
vasp_results_dir = '../Database/TestDB/second_run/VASP6/BEXT=0.001/outdir'

screener_after(wdatalist)


