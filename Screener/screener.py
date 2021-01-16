import pandas as pd
import os
import numpy as np
import tqdm
import tarfile
import math
import matplotlib.pyplot as plt
from pymatgen.io.vasp import Poscar
from pymatgen.symmetry.analyzer import SpacegroupAnalyzer
### LOCAL IMPORTS ###
from Tools import POSCAR_reader
import time
import timeit
import scipy
from collections import Counter


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
    df['polymorphs'] = 1 + (df.groupby('compound').compound.transform(lambda x: x.duplicated(keep=False).sum()))  # check how many entries still have the same composition - are polymorphs with diffeent structures - and count them
    df.to_csv(datalist.replace('.csv', '_no.duplicates.csv'))


def mag_sites_calculator(ID):
    """Determines how many unique magnetic sites are present in the structure.
    Takes ID as input, looks up structure file in the datadir and returns number of unique sites as integer"""

    ### List of Magnetic Atoms
    magnetic = ['Ti', 'V', 'Cr', 'Mn', 'Fe', 'Co', 'Ni',  # Cu, Sc, Y, and everything else is considered not magnetic
                 'Nb', 'Mo', 'Ru', 'Rh', 'Pd',  # we sieve out Cd at earlier step, so I could have omitted it from this list but this approach is more "universal"
                 'Ce', 'Pr', 'Nd', 'Sm', 'Eu', 'Gd', 'Tb', 'Dy', 'Ho', 'Er', 'Tm', 'Yb']
    bad_notation = False                                                    # sometimes aflow files have alphabet instead of element symbols in structure POSCAR I call it bad notation
    rename_list = ['A', 'B', 'C', 'D', 'E', 'F', 'G', 'H', 'I', 'J', 'K']   # A list for fixing bad notation
                                                                            # I don't expect there to be a composition with more than 11 distinctive components so here I stop at K

    eltable = np.empty([0, 2])
    df = pd.read_csv(wdatalist, index_col=0, sep=',')

    structure_file = wdatadir_structure + str(ID)         # still need to decide weather I should

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
                        if bad_notation == True: # Apply correction for bad notation replacing all alphabet letters with proper chemical symbols taken from list of elements in composition
                            n = rename_list.index(elem)
                            temp_line = (df.loc[ID, 'species']).strip('[').strip(']').strip(' ')
                            temp_line2 = temp_line.replace("'", "")
                            species_list = np.asarray(temp_line2.split('; '))
                            elem = species_list[n]
                        eltable = np.append(eltable, [[elem, wyckoff]], axis=0)
    # This section with masks probably could be done better and shorter but with such small arrays it shouldn't matter much
    mask = np.in1d(eltable[:, 0], magnetic)  # Check what entries in Wyckoff list are in list of magnetic atoms
    eltable2 = eltable[mask]                 # new array of only magnetic ones using previous mask
    unique_keys, mask2 = np.unique(eltable2[:, 1], return_index=True)  # check what entries of the new array have unique wyckoff sites - creates a list of indicies corresponding to said entries
    site_counter = np.size(eltable2[mask2], axis=0)                    # get number of unique sites counting number of obtained indicies
    return site_counter


def mag_sites_calculator_MP(ID):
    """Determines how many unique magnetic sites are present in the structure.
    Takes ID as input, looks up structure file in the datadir and returns number of unique sites as integer"""

    ### List of Magnetic Atoms
    magnetic = ['Ti', 'V', 'Cr', 'Mn', 'Fe', 'Co', 'Ni',        # Cu, Sc, Y, and everything else is considered not magnetic
                 'Nb', 'Mo', 'Ru', 'Rh', 'Pd',                  # we sieve out Cd at earlier step, so I could have omitted it from this list but this approach is more "universal"
                 'Ce', 'Pr', 'Nd', 'Sm', 'Eu', 'Gd', 'Tb', 'Dy', 'Ho', 'Er', 'Tm', 'Yb']

    eltable = np.empty([0, 2])
    structure_file = wdatadir_structure + str(ID) + '.cif'
    poscar_content = POSCAR_reader.read(structure_file)
    poscar_string = ''.join(poscar_content)       # merging poscar content in a single string so pymatgen can read it
    poscar = Poscar.from_string(poscar_string)    # using pymatgen to acquire our structure from poscar content
    structure = poscar.structure
    analyzed_structure = str(SpacegroupAnalyzer(structure).get_symmetrized_structure()).split('\n')
    for num,  line in enumerate(analyzed_structure):
        if num > 7:
            elem = line.split()[1]  # element symbol is after Wyckoff x y z coordinates so 4th item in list
            wyckoff = line.split()[5]
            eltable = np.append(eltable, [[elem, wyckoff]], axis=0)
    mask = np.in1d(eltable[:, 0], magnetic)  # Check what entries in Wyckoff list are in list of magnetic atoms
    eltable2 = eltable[mask]                 # new array of only magnetic ones using previous mask
    unique_keys, mask2 = np.unique(eltable2[:, 1], return_index=True)  # check what entries of the new array have unique wyckoff sites - creates a list of indicies corresponding to said entries
    site_counter = np.size(eltable2[mask2], axis=0)                    # get number of unique sites counting number of obtained indicies
    return site_counter


def composition_analysis(ID):
    """ """
    rare_earths = ['Sc', 'Y', 'La', 'Ce', 'Pr', 'Nd', 'Pm', 'Sm', 'Eu', 'Gd', 'Tb', 'Dy', 'Ho', 'Er', 'Tm', 'Yb', 'Lu']
    oxygen = ['O']
    rare_total = 0
    oxygen_total = 0

    structure_file = wdatadir_structure + str(ID) + '.cif'
    poscar_content = POSCAR_reader.read(structure_file)
    poscar_string = ''.join(poscar_content)       # merging poscar content in a single string so pymatgen can read it
    poscar = Poscar.from_string(poscar_string)    # using pymatgen to acquire our structure from poscar content
    structure = poscar.structure

    number = len(structure.cart_coords)

    el_list = [str(i).replace('Element ', '') for i in structure.species]

    c = Counter(el_list)

    for i in c.keys():
        if i in oxygen:
            oxygen_total = oxygen_total + c[i]
    o_content = 100*oxygen_total/number

    for i in c.keys():
        if i in rare_earths:
            rare_total = rare_total + c[i]
    rare_content = 100*rare_total/number

    return number, rare_content, o_content


def timing(ID):
    """Reads OUTCARs and returns time taken to calculate entry in minutes (total: undeformed + all deformations)"""

    time_spent = 0
    deformations = os.listdir(vasp_results_dir + '/' + str(ID))
    for deformation in deformations:
        with open(vasp_results_dir + '/' + str(ID) + '/' + deformation + "/OUTCAR") as search:
            for line in search:
                if 'User time' in line:
                    time_spent = time_spent + (float(line.split()[-1]))
    return round(time_spent/60)


def completion_date(ID):
    """Reads OUTCARs and returns the date all calculations for an entry were finished
    (latest among undeformed and all deformations)"""

    # strictly speaking using time module here is an overkill.
    # Just sorting strings should be sufficient as VASP gives them in suitable format already.
    # this is more flexible and timing for both approaches are similar so...

    dates = []
    deformations = os.listdir(vasp_results_dir + '/' + str(ID))
    for deformation in deformations:
        with open(vasp_results_dir + '/' + str(ID) + '/' + deformation + "/OUTCAR") as search:
            for line in search:
                if 'executed on' in line:
                    dates.append(time.strptime((line.split()[-2]), "%Y.%m.%d"))
    date = time.strftime('%Y-%m-%d', max(dates))
    return date


def memory_used(ID):
    """Reads OUTCARs and returns memmory used to calculate entry (total: undeformed + all deformations)"""

    memory = 0
    deformations = os.listdir(vasp_results_dir + '/' + str(ID))
    for deformation in deformations:
        with open(vasp_results_dir + '/' + str(ID) + '/' + deformation + "/OUTCAR") as search:
            for line in search:
                if 'Maximum memory used' in line:
                    memory = memory + (float(line.split()[-1]))
    memory = memory * 9.5367432e-7  # converting kb to Gb
    return round(memory, 2)


def geometry_after(ID, deformation):
    """Reads OUTCARs and returns lattice parameters as detrmined by vasp.
    If chosen deformation is n/a for this ID an empty string will be returned.
    WARNING! The parameters returned are for whaterver supercell vasp was working with.
    They may not be the actual lattice parameters for the compound but instead multiples of said parameters.
    Use with caution."""

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
        return '', '', ''


def moment_volume_after(ID='', deformation='', file=''):
    """Read OUTCAR file for chosen deformation type and return moment, volume and
    value of magnetic field based on moment and volume.
    If chosen deformation is n/a for this ID an empty string will be returned
    If file is specified all information will be taken from it instead"""

    if file == '':
        path_to_data = vasp_results_dir + '/' + str(ID) + '/' + deformation + "/OUTCAR"

    else:
        path_to_data = file

    with open(path_to_data) as search:
        moment_lines = []
        for line in search:
            if 'volume of cell' in line:
                volume = (float(line.split()[-1]))
            if 'tot' in line:
                moment_lines.append(line)
    moment = float(moment_lines[-2].split()[-1])
    field_after = calculate_mag_field(moment, volume)
    return moment, volume, field_after



def sym_after(ID, deformation):
    """ Gives symmetry group determined by VASP for chosen deformation. If Deformation is n/a returns an empty string"""
    if deformation in os.listdir(vasp_results_dir + '/' + str(ID)):
        with open(vasp_results_dir + '/' + str(ID) + '/' + deformation + "/OUTCAR") as search:
            for line in search:
                if 'Routine SETGRP: Setting up the symmetry group for a' in line:
                    symmetry = (next(search)).strip("\n")
        return symmetry
    else:
        return ''


def energy(ID, deformation):
    """ Gives total energy calculated by VASP for chosen deformation. If Deformation is n/a returns an empty string"""
    if deformation in os.listdir(vasp_results_dir + '/' + str(ID)):
        with open(vasp_results_dir + '/' + str(ID) + '/' + deformation + "/OUTCAR") as search:
            energy_lines = []
            for line in search:
                if 'free  energy   TOTEN' in line:
                    energy_lines.append(line)
        tot_energy = round(float(energy_lines[-1].split()[-2]), 4)
        return tot_energy
    else:
        return ''


def archive_reader(ID, deformation):
    """Function to work with archived files created from rundir
    At the moment reads number of RMM steps from OSZICAR.
    To be repuporsed later for other queries"""

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
        return ''


def mag_sites_difference(ID):
    # Works after calculation, reads OUTCAR and checks values of moments for different magnetic sites and returns biggest
    # difference between moments of sublattices or just all significantly different moments...
    return


def linfit(x, y):
    coeffs = np.polyfit(x, y, 1)
    poly1d_fn = np.poly1d(coeffs)

    # fit quality
    yhat = poly1d_fn(x)  # or [p(z) for z in x]
    ybar = np.sum(y) / len(y)  # or sum(y)/len(y)
    ssreg = np.sum((yhat - ybar) ** 2)  # or sum([ (yihat - ybar)**2 for yihat in yhat])
    sstot = np.sum((y - ybar) ** 2)  # or sum([ (yi - ybar)**2 for yi in y])
    R_square = ssreg / sstot
    return coeffs[0], R_square


def parabfit(x, y):
    coeffs = np.polyfit(x, y, 2)
    poly1d_fn = np.poly1d(coeffs)
    crit = poly1d_fn.deriv().r
    extremum_x = float(crit[crit.imag == 0].real)
    extremum_y = poly1d_fn(extremum_x)
    extremum = [extremum_x, extremum_y]
    # test = poly1d_fn.deriv(2)(r_crit) # can check if min or max if necessary
    return extremum


def magnetoelastic(ID, deformation):
    if deformation + '_inc' in os.listdir(vasp_results_dir + '/' + str(ID)):
        moment_inc, volume_inc, field_inc = moment_volume_after(ID, deformation+'_inc')
        moment_dec, volume_dec, field_dec = moment_volume_after(ID, deformation + '_dec')
        moment, volume, field = moment_volume_after(ID, 'undeformed')
        y = [field_dec, field, field_inc]
        x = [0.95, 1, 1.05]

        slope, fit_quality = linfit(x, y)

        if fit_quality < 0.75:
            extremum = parabfit(x, y)

            if abs(extremum[0] - 0.95) > abs(1.05 - extremum[0]):
                x2 = [0.95, 1, extremum[0]]
                y2 = [field_dec, field, extremum[1]]
                new_slope, new_fit_quality = linfit(x2, y2)
            else:
                x2 = [extremum[0], 1, 1.05]
                y2 = [extremum[1], field, field_inc]
                new_slope, new_fit_quality = linfit(x2, y2)
        else:
            new_slope, new_fit_quality = slope, ''

        return slope, fit_quality, new_slope, new_fit_quality
    else:
        return '', '', '', ''


def magnetoelastic_moment(moment_dec, moment, moment_inc):

    if moment and moment_dec and moment_inc != 0.0:
        moment_change_dec = moment_dec/moment
        moment_change_inc = moment_inc/moment
    else:
        return 0, 0, 0, 0

    y = [moment_change_dec, 1, moment_change_inc]
    x = [0.95, 1, 1.05]

    slope, fit_quality = linfit(x, y)

    if fit_quality < 0.75:
        extremum = parabfit(x, y)

        if abs(extremum[0] - 0.95) > abs(1.05 - extremum[0]):
            x2 = [0.95, 1, extremum[0]]
            y2 = [moment_change_dec, 1, extremum[1]]
            new_slope, new_fit_quality = linfit(x2, y2)
        else:
            x2 = [extremum[0], 1, 1.05]
            y2 = [extremum[1], 1, moment_change_inc]
            new_slope, new_fit_quality = linfit(x2, y2)
    else:
        new_slope, new_fit_quality = slope, ''

    return slope, fit_quality, new_slope, new_fit_quality


def screener_before(datalist, database_type):
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

            number_atoms, rare_content, oxygen_content = composition_analysis(item)

        # Write number of atoms into datalist
            df.loc[item, 'number_atoms'] = number_atoms
        # Write rare earth content into datalist
            df.loc[item, 'number_atoms'] = rare_content
        #   Write oxygen content into datalist
            df.loc[item, 'oxygen_content'] = oxygen_content

        # Write number of sites into datalist:
            if database_type== 'aflow':
                df.loc[item, 'mag_sites'] = mag_sites_calculator(item)
            elif database_type== 'MP':
                df.loc[item, 'mag_sites'] = mag_sites_calculator_MP(item)

        # Write internal magnetic field into datalist:
            df.loc[item, 'mag_field'] = calculate_mag_field(df.loc[item, 'moment_cell'], df.loc[item, 'volume_cell'])

        # Wrire an updated datalist back to file
        df.to_csv(datalist.replace('.csv', '_updated.csv'))

        # drops everything that does not fit the criteria, creating new datalist file at each sieve
    sieve(datalist.replace('.csv', '_updated.csv'), 'mag_field', min_mag_field)
    sieve(datalist.replace('.csv', '_updated_sieved.mag.field.csv'), 'mag_sites', min_site_number)

        # Drops duplicates, creates new datalist file __ PERHAPS THIS SHOULD BE FIRST STEP TO LOWER FURTHER WORKLOAD
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


def moment_volume_after_temp(ID, deformation):
    """Read OUTCAR file for chosen deformation type and return moment, volume and
    value of magnetic field based on moment and volume.
    If chosen deformation is n/a for this ID an empty string will be returned"""

    if deformation in os.listdir(vasp_results_BEXT + '/' + str(ID)):
        with open(vasp_results_BEXT + '/' + str(ID) + '/' + deformation + "/OUTCAR") as search:
            moment_lines = []
            for line in search:
                if 'volume of cell' in line:
                    volume = (float(line.split()[-1]))
                if 'tot' in line:
                    moment_lines.append(line)
        moment = float(moment_lines[-2].split()[-1])
        field_after = calculate_mag_field(moment, volume)
        return moment, volume, field_after
    else:
        return '', '', ''


def do_datalist(datalist):
    """A small instance of screener, just to apply some parameter to the chosen datalist"""

    df = pd.read_csv(datalist, index_col=0, sep=',')
    with tqdm.tqdm(total=len(df.index)) as pbar:        # A wrapper that creates nice progress bar
        pbar.set_description("Processing datalist")
        for item in df.index.tolist():
            pbar.update(1)                              # Updating progress bar at each step
            # Write number of atoms into datalist
            a = df.loc[item, 'Mel_a_fit']
            b = df.loc[item, 'Mel_b_fit']
            c = df.loc[item, 'Mel_c_fit']

            if a == '':
                a = 1
            if b == '':
                b = 1
            if c == '':
                c = 1

            poorest_fit = min(a, b, c)
            df.loc[item, 'poorst_fit'] = poorest_fit

    df.to_csv(datalist.replace(".csv", '_new' + '.csv'))


def sieve_temp01(datalist):
    """Removes entries from datalist according to selected sieving criteria"""
    df = pd.read_csv(datalist, index_col=0, sep=',')
    with tqdm.tqdm(total=len(df.index)) as pbar:        # A wrapper that creates nice progress bar
        pbar.set_description("Processing datalist")

        typelist = ['a_inc_fail_reason',
                    'undeformed_fail_reason',
                    'c_dec_fail_reason',
                    'b_dec_fail_reason',
                    'a_dec_fail_reason',
                    'c_inc_fail_reason',
                    'b_inc_fail_reason']

        line = 'run failed, set accuracy not reached: dE is'

        for item in df.index.tolist():
            pbar.update(1)                              # Updating progress bar at each step
            e_list = []

            for element in typelist:
                reason = str(df.loc[item, element])
                if line in reason:

                    dEi = reason.split('is ')[1].split(',')[0]
                    if dEi != '':
                        # print(dEi)
                        e_list.append(float(dEi))

            if type(e_list) != type(None):
                df.loc[item, 'dE'] = max(e_list)

    # df2.to_csv(datalist.replace(".csv", '_O.csv'))
    df.to_csv(datalist.replace(".csv", "_dE.csv"))


def sieve_temp02(datalist):
    """Removes entries from datalist according to selected sieving criteria"""

    df = pd.read_csv(datalist, index_col=0, sep=',')
    df2 = pd.DataFrame(columns=df.columns)

    typelist = ['a_inc_fail_reason',
                'undeformed_fail_reason',
                'c_dec_fail_reason',
                'b_dec_fail_reason',
                'a_dec_fail_reason',
                'c_inc_fail_reason',
                'b_inc_fail_reason']

    line = 'run failed, set accuracy not reached: dE is'
    for item in df.index.tolist():
        for element in typelist:
            if line in str(df.loc[item, element]):
                rows = df.loc[item]
                df2 = df2.append(rows)
                df = df.drop([item], axis=0)
                break


    df2 = df2.drop_duplicates()
    df2.to_csv(datalist.replace(".csv", '_accuracy.csv'))
    df.to_csv(datalist.replace(".csv", "_except_accuracy.csv"))


def OUTCAR_reader(ID='', folder=''):
    """A universdal function to get the whole content of the OUTCAR file for certain ID
    Speeds up the screening process, but lacks flexibility of separate functions
    If folder is specified will look at all deformations for that single folder"""

    if folder == '':
        path_to_data = vasp_results_dir + '/' + str(ID)
    else:
        path_to_data = folder

    moment_lines = []
    energy_lines = []

    def_list = ['undeformed', 'a_dec', 'a_inc', 'b_dec', 'b_inc', 'c_dec', 'c_inc']

    deformations = os.listdir(path_to_data)

    dd = pd.DataFrame(0, index=def_list, dtype='float', columns=['time_spent', 'date_complete', 'memory', 'volume', 'moment', 'symmetry', 'tot_energy'])
    dd[['date_complete', 'symmetry']] = dd[['date_complete', 'symmetry']].astype('str')

    for deformation in deformations:
        with open(path_to_data + '/' + deformation + "/OUTCAR") as search:
            for line in search:

                if 'User time' in line:
                    dd.at[deformation, 'time_spent'] = (float(line.split()[-1]))/60

                if 'executed on' in line:
                    dd.loc[deformation, 'date_complete'] = time.strftime('%Y-%m-%d', (time.strptime((line.split()[-2]), "%Y.%m.%d")))

                if 'Maximum memory used' in line:
                    dd.at[deformation, 'memory'] = (float(line.split()[-1]))

                if 'volume of cell' in line:
                    dd.at[deformation, 'volume'] = (float(line.split()[-1]))

                if 'Routine SETGRP: Setting up the symmetry group for a' in line:
                    dd.loc[deformation, 'symmetry'] = (next(search)).strip("\n").split()[1]

                if 'tot' in line:
                    moment_lines.append(line)

                if 'free  energy   TOTEN' in line:
                    energy_lines.append(line)

        dd.at[deformation, 'moment'] = float(moment_lines[-2].split()[-1])
        dd.at[deformation, 'tot_energy'] = round(float(energy_lines[-1].split()[-2]), 4)


    # time_spent = round(time_spent / 60)
    # date = time.strftime('%Y-%m-%d', max(dates))
    # memory = round(memory * 9.5367432e-7, 2)  # converting kb to Gb and rounding
    return dd


def screener_after(datalist):
    """Second main function used for applying criteria obtained after calculations and performing corresponding checks"""
    df = pd.read_csv(datalist, index_col=0, sep=',')

    columns_to_remove = ['a_inc', 'a_dec', 'b_inc', 'b_dec', 'c_inc', 'c_dec', 'undeformed',
                         'a_inc_fail_reason', 'a_dec_fail_reason', 'b_inc_fail_reason', 'b_dec_fail_reason', 'c_dec_fail_reason', 'c_inc_fail_reason', 'undeformed_fail_reason']
    for i in columns_to_remove:
        try:
            df.drop(columns_to_remove, inplace=True, axis=1)
        except:
            pass

    with tqdm.tqdm(total=len(df.index)) as pbar:        # A wrapper that creates nice progress bar
        pbar.set_description("Processing datalist")
        for item in df.index.tolist():
            pbar.update(1)                              # Updating progress bar at each step
            OUTCAR_results = OUTCAR_reader(item)

            number_atoms, rare_content, oxygen_content = composition_analysis(item)

            df.at[item, 'number_atoms'] = number_atoms
            # Rare earth content into datalist
            df.at[item, 'rare_content'] = rare_content
            # Oxygen content into datalist
            df.at[item, 'oxygen_content'] = oxygen_content

            # Memory used for calculation
            df.at[item, 'memory_used'] = round(OUTCAR_results['memory'].sum() * 9.5367432e-7, 2)
            # Time used for calculation
            df.at[item, 'time_to_calculate'] = round(OUTCAR_results['time_spent'].sum() / 60, 2)
            # Date calculation was finished on
            df.at[item, 'date_complete'] = OUTCAR_results['date_complete'].sort_values()[-1]
            # Total Energy
            df.at[item, 'energy'] = OUTCAR_results.at['undeformed', 'tot_energy']

            # Volume after
            df.at[item, 'moment_u'] = OUTCAR_results.at['undeformed', 'moment']
            # Moment after
            df.at[item, 'volume_u'] = OUTCAR_results.at['undeformed', 'volume']
            # Mag_field after
            df.at[item, 'magF_u'] = abs(calculate_mag_field(OUTCAR_results.at['undeformed', 'moment'], OUTCAR_results.at['undeformed', 'volume']))

            # Magneto Elastic parameter for all axes
            df.at[item, 'Mel_a'] = magnetoelastic_moment(OUTCAR_results.at['a_dec', 'moment'], OUTCAR_results.at['undeformed', 'moment'], OUTCAR_results.at['a_inc', 'moment'])[2]
            df.at[item, 'Mel_b'] = magnetoelastic_moment(OUTCAR_results.at['b_dec', 'moment'], OUTCAR_results.at['undeformed', 'moment'], OUTCAR_results.at['b_inc', 'moment'])[2]
            df.at[item, 'Mel_c'] = magnetoelastic_moment(OUTCAR_results.at['c_dec', 'moment'], OUTCAR_results.at['undeformed', 'moment'], OUTCAR_results.at['c_inc', 'moment'])[2]

            df.at[item, 'Mel_V'] = magnetoelastic_moment(OUTCAR_results.at['V_dec', 'moment'], OUTCAR_results.at['undeformed', 'moment'], OUTCAR_results.at['V_inc', 'moment'])[2]

        # Total Magneto Elastic parameter
        df['Mel_full'] = (df['Mel_a'] ** 2 + df['Mel_b'] ** 2 + df['Mel_c'] ** 2) ** 0.5

    df.to_csv(datalist.replace(".csv", '_out' + '.csv'))

# ------------------------------------------------------------------------------------------------------- #
#  _____                                           _       _____ _             _     _   _                #
# /  __ \                                         | |     /  ___| |           | |   | | | |               #
# | /  \/ ___  _ __ ___  _ __ ___   __ _ _ __   __| |___  \ `--.| |_ __ _ _ __| |_  | |_| | ___ _ __ ___  #
# | |    / _ \| '_ ` _ \| '_ ` _ \ / _` | '_ \ / _` / __|  `--. \ __/ _` | '__| __| |  _  |/ _ \ '__/ _ \ #
# | \__/\ (_) | | | | | | | | | | | (_| | | | | (_| \__ \ /\__/ / || (_| | |  | |_  | | | |  __/ | |  __/ #
#  \____/\___/|_| |_| |_|_| |_| |_|\__,_|_| |_|\__,_|___/ \____/ \__\__,_|_|   \__| \_| |_/\___|_|  \___| #
#                                                                                                         #
# ------------------------------------------------------------------------------------------------------- #

### Setting which database we work with ###
# wdatadir_structure = '../Database/aflow/datadir_structure_relaxed/'
# wdatalist = '../Database/aflow/datalist_updated_sieved.mag.field_sieved.mag.sites.csv'


wdatadir_structure = 'D:/MCES/MP/datadir/'
wdatalist = 'D:/MCES/MP/step4_success_sieved.csv'

outdir = 'D:/MCES/MP/outdir_VVV'
vasp_results_dir = outdir
vasp_results_BEXT = 'D:/MCES/MP/outdir_BEXT_1/out'
# screener_before('X:/MCES/MP/step2.csv', 'MP')


# screener_after('step4_success_sieved.csv')
screener_after('D:/MCES/MP/V_beforeRun_afterRun_success_sieved.csv')

# desired_width = 320
# pd.set_option("display.max_rows", None, "display.max_columns", None)
# pd.set_option('display.width', desired_width)
# pd.set_option('display.max_columns', 10)

# print(moment_volume_after('mp-778', 'a_inc'))
# something is wrong with this one, also I don't trust deformations again are they done rightly????