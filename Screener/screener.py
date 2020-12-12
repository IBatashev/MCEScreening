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


def mag_sites_calculator(ID):
    """Determines how many unique magnetic sites are present in the structure.
    Takes ID as input, looks up structure file in the datadir and returns number of unique sites as integer"""

    ### List of Magnetic Atoms
    magnetic = [ 'Ti', 'V', 'Cr', 'Mn', 'Fe', 'Co', 'Ni',  # Cu, Sc, Y, and everythong else is considered not magnetic
                 'Nb', 'Mo', 'Ru', 'Rh', 'Pd', # we sieve out Cd at earliear step, so I could have ommitted it from this list but this approach is more "universal"
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
                        if bad_notation == True:# Apply correction for bad notation replacing all alphabet letters with proper chemical symbols taken from list of elements in composition
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
    magnetic = [ 'Ti', 'V', 'Cr', 'Mn', 'Fe', 'Co', 'Ni', #'Cu',  # Sc, Y, and everythong else is considered not magnetic
                 'Nb', 'Mo', 'Ru', 'Rh', 'Pd', # we sieve out Cd at earliear step, so I could have ommitted it from this list but this approach is more "universal"
                 'Ce', 'Pr', 'Nd', 'Sm', 'Eu', 'Gd', 'Tb', 'Dy', 'Ho', 'Er', 'Tm', 'Yb']

    eltable = np.empty([0, 2])
    # df = pd.read_csv(wdatalist, index_col=0, sep=',')
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


def number_of_atoms(ID):
    """counts how many atoms are listed in poscar for chosen structure"""
    structure_file = wdatadir_structure + str(ID) + '.cif'
    poscar_content = POSCAR_reader.read(structure_file)
    poscar_string = ''.join(poscar_content)       # merging poscar content in a single string so pymatgen can read it
    poscar = Poscar.from_string(poscar_string)    # using pymatgen to acquire our structure from poscar content
    structure = poscar.structure
    number = len(structure.cart_coords)
    return number


def rare_earth_content(ID):
    """Calculates % of rare earth atoms from total atom number as well as total number of atoms"""

    rare_earths = ['Sc', 'Y', 'La', 'Ce', 'Pr', 'Nd', 'Pm', 'Sm', 'Eu', 'Gd', 'Tb', 'Dy', 'Ho', 'Er', 'Tm', 'Yb', 'Lu']

    structure_file = wdatadir_structure + str(ID) + '.cif'
    poscar_content = POSCAR_reader.read(structure_file)
    poscar_string = ''.join(poscar_content)       # merging poscar content in a single string so pymatgen can read it
    poscar = Poscar.from_string(poscar_string)    # using pymatgen to acquire our structure from poscar content
    structure = poscar.structure

    number = len(structure.cart_coords)
    el_list = [str(i).replace('Element ', '') for i in structure.species]
    c = Counter(el_list)
    rare_total = 0
    for i in c.keys():
        if i in rare_earths:
            rare_total = rare_total + c[i]
    rare_content = 100*rare_total/number

    return rare_content, number


def oxygen_content(ID):
    """Calculates % of oxygen atoms from total atom number as well as total number of atoms"""

    oxygen = ['O']

    structure_file = wdatadir_structure + str(ID) + '.cif'
    poscar_content = POSCAR_reader.read(structure_file)
    poscar_string = ''.join(poscar_content)       # merging poscar content in a single string so pymatgen can read it
    poscar = Poscar.from_string(poscar_string)    # using pymatgen to acquire our structure from poscar content
    structure = poscar.structure

    number = len(structure.cart_coords)
    el_list = [str(i).replace('Element ', '') for i in structure.species]
    c = Counter(el_list)
    oxygen_total = 0
    for i in c.keys():
        if i in oxygen:
            oxygen_total = oxygen_total + c[i]
    o_content = 100*oxygen_total/number

    return o_content


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


def timing(ID):
    """Reads OUTCARs and returns time taken to calculate entry in minutes (total: undeformed + all deformations)"""

    time = 0
    deformations = os.listdir(vasp_results_dir + '/' + str(ID))
    for deformation in deformations:
        with open(vasp_results_dir + '/' + str(ID) + '/' + deformation + "/OUTCAR") as search:
            for line in search:
                if 'User time' in line:
                    time = time + (float(line.split()[-1]))
    return round(time/60)


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


def moment_volume_after(ID, deformation):
    """Read OUTCAR file for chosen deformation type and return moment, volume and
    value of magnetic field based on moment and volume.
    If chosen deformation is n/a for this ID an empty string will be returned"""

    if deformation in os.listdir(vasp_results_dir + '/' + str(ID)):
        with open(vasp_results_dir + '/' + str(ID) + '/' + deformation + "/OUTCAR") as search:
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


def mag_sites_difference(ID):
    # Works after calculation, reads OUTCAR and checks values of moments for different magnetic sites and returns biggest
    # difference between moments of sublattices or just all significantly different moments...
    return


# def magnetoelastic(ID, deformation, plot=False):
#     """"""
#
#     if deformation + '_inc' in os.listdir(vasp_results_dir + '/' + str(ID)):
#         moment_inc, volume_inc, field_inc = moment_volume_after(ID, deformation+'_inc')
#         moment_dec, volume_dec, field_dec = moment_volume_after(ID, deformation + '_dec')
#         moment, volume, field = moment_volume_after(ID, 'undeformed')
#         y = [field_dec, field, field_inc]
#         x = [0.95, 1, 1.05]
#         coeff = np.polyfit(x, y, 1)
#         if plot == True:
#             print('Plotting...')
#             poly1d_fn = np.poly1d(coeff)
#             plt.plot(x, y, 'yo', x, poly1d_fn(x), '--k')
#             plt.xlim(0.9, 1.1)
#             plt.show()
#         return coeff[0] #, field_dec, field_inc
#     else:
#         return '' #, '', ''

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


def magnetoelastic_moment(ID, deformation):
    if deformation + '_inc' in os.listdir(vasp_results_dir + '/' + str(ID)):
        moment_inc, volume_inc, field_inc = moment_volume_after(ID, deformation+'_inc')
        moment_dec, volume_dec, field_dec = moment_volume_after(ID, deformation + '_dec')
        moment, volume, field = moment_volume_after(ID, 'undeformed')
        y = [moment_dec/moment, 1, moment_inc/moment]
        x = [0.95, 1, 1.05]

        slope, fit_quality = linfit(x, y)

        if fit_quality < 0.75:
            extremum = parabfit(x, y)

            if abs(extremum[0] - 0.95) > abs(1.05 - extremum[0]):
                x2 = [0.95, 1, extremum[0]]
                y2 = [moment_dec/moment, 1, extremum[1]]
                new_slope, new_fit_quality = linfit(x2, y2)
            else:
                x2 = [extremum[0], 1, 1.05]
                y2 = [extremum[1], 1, moment_inc/moment]
                new_slope, new_fit_quality = linfit(x2, y2)
        else:
            new_slope, new_fit_quality = slope, ''

        return slope, fit_quality, new_slope, new_fit_quality
    else:
        return '', '', '', ''


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

            rare_content, number_atoms = rare_earth_content(item)

        # Write number of atoms into datalist
            df.loc[item, 'number_atoms'] = number_atoms
        # Write rare earth content into datalist
            df.loc[item, 'number_atoms'] = rare_content
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


def screener_after(datalist):
    """Second main function used for applying criteria obtained after calculations and performing corresponding checks"""

    min_mag_field = 0.45

    df = pd.read_csv(datalist, index_col=0, sep=',')
    with tqdm.tqdm(total=len(df.index)) as pbar:        # A wrapper that creates nice progress bar
        pbar.set_description("Processing datalist")
        for item in df.index.tolist():
            pbar.update(1)                              # Updating progress bar at each step

            rare_content, number_atoms = rare_earth_content(item)

            # Write number of atoms into datalist
            df.loc[item, 'number_atoms'] = number_atoms

            # Write rare earth content into datalist
            df.loc[item, 'rare_content'] = rare_content

            # Write oxygen content into datalist
            df.loc[item, 'oxygen_content'] = oxygen_content(item)

            # Get memory used for calculation
            df.loc[item, 'memory_used'] = memory_used(item)

            # Get time used for calculation
            df.loc[item, 'time_to_calculate'] = timing(item)

            # Get date calculation was finished on
            df.loc[item, 'date_complete'] = completion_date(item)

            moment_u, volume_u, magF_u = moment_volume_after(item, 'undeformed')
            # Get volume after
            # df.loc[item, 'moment_u'] = moment_volume_after(item, 'undeformed')[0]
            df.loc[item, 'moment_u'] = moment_u

            # Get moment after
            # df.loc[item, 'volume_u'] = moment_volume_after(item, 'undeformed')[1]
            df.loc[item, 'volume_u'] = volume_u

            # Calculate mag_field after
            # df.loc[item, 'magF_u'] = abs(moment_volume_after(item, 'undeformed')[2]) # we are only interested in absolute value
            df.loc[item, 'magF_u'] = abs(magF_u)

            # Get primitive cell parameters determined by vasp for this structure
            df.loc[item, 'VASP_a_param'] = geometry_after(item, 'undeformed')[0]
            df.loc[item, 'VASP_b_param'] = geometry_after(item, 'undeformed')[1]
            df.loc[item, 'VASP_c_param'] = geometry_after(item, 'undeformed')[2]

            # Get algorithm step at which optimization stopped
            df.loc[item, 'step'] = archive_reader(item, 'undeformed')

            # Get energy
            df.loc[item, 'energy'] = energy(item, 'undeformed')

            # Calculate magnetoelastic parameter initial with fit quality and adjusted with fit quality
            Mel_a, Mel_a_fit, Mel_aa, Mel_aa_fit = magnetoelastic(item, 'a')
            Mel_b, Mel_b_fit, Mel_bb, Mel_bb_fit = magnetoelastic(item, 'b')
            Mel_c, Mel_c_fit, Mel_cc, Mel_cc_fit = magnetoelastic(item, 'c')

            # Calculate magnetoelastic parameter via moments initial with fit quality and adjusted with fit quality
            Mel_mom_a, Mel_mom_a_fit, Mel_mom_aa, Mel_mom_aa_fit = magnetoelastic_moment(item, 'a')
            Mel_mom_b, Mel_mom_b_fit, Mel_mom_bb, Mel_mom_bb_fit = magnetoelastic_moment(item, 'b')
            Mel_mom_c, Mel_mom_c_fit, Mel_mom_cc, Mel_mom_cc_fit = magnetoelastic_moment(item, 'c')


            # OPTIONAL TEMPORARY FOR EXTRA ANALYSIS
            df.loc[item, 'a_inc_m'] = moment_volume_after(item, 'a_inc')[0]
            df.loc[item, 'a_dec_m'] = moment_volume_after(item, 'a_dec')[0]
            df.loc[item, 'b_inc_m'] = moment_volume_after(item, 'b_inc')[0]
            df.loc[item, 'b_dec_m'] = moment_volume_after(item, 'b_dec')[0]
            df.loc[item, 'c_inc_m'] = moment_volume_after(item, 'c_inc')[0]
            df.loc[item, 'c_dec_m'] = moment_volume_after(item, 'c_dec')[0]

            df.loc[item, 'Mel_a'] = Mel_aa
            df.loc[item, 'Mel_b'] = Mel_bb
            df.loc[item, 'Mel_c'] = Mel_cc
            df.loc[item, 'Mel_a_fit'] = Mel_a_fit
            df.loc[item, 'Mel_b_fit'] = Mel_b_fit
            df.loc[item, 'Mel_c_fit'] = Mel_c_fit

            df.loc[item, 'Mel_mom_a'] = Mel_mom_aa
            df.loc[item, 'Mel_mom_b'] = Mel_mom_bb
            df.loc[item, 'Mel_mom_c'] = Mel_mom_cc
            df.loc[item, 'Mel_mom_a_fit'] = Mel_mom_a_fit
            df.loc[item, 'Mel_mom_b_fit'] = Mel_mom_b_fit
            df.loc[item, 'Mel_mom_c_fit'] = Mel_mom_c_fit

            if Mel_a == '':
                Mel_a = 0
                Mel_aa = 0
            if Mel_b == '':
                Mel_b = 0
                Mel_bb = 0
            if Mel_c == '':
                Mel_c = 0
                Mel_cc = 0

            df.loc[item, 'Mel_full_old'] = math.sqrt((Mel_a**2) + (Mel_b**2) + (Mel_c**2))
            df.loc[item, 'Mel_full'] = math.sqrt((Mel_aa**2) + (Mel_bb**2) + (Mel_cc**2))


            if Mel_mom_a == '':
                Mel_mom_a = 0
                Mel_mom_aa = 0
            if Mel_mom_b == '':
                Mel_mom_b = 0
                Mel_mom_bb = 0
            if Mel_mom_c == '':
                Mel_mom_c = 0
                Mel_mom_cc = 0

            df.loc[item, 'Mel_mom_full_old'] = math.sqrt((Mel_mom_a**2) + (Mel_mom_b**2) + (Mel_mom_c**2))
            df.loc[item, 'Mel_mom_full'] = math.sqrt((Mel_mom_aa**2) + (Mel_mom_bb**2) + (Mel_mom_cc**2))

            try:
                magF_a = (moment_volume_after_temp(item, 'Applied_Field')[2])
                df.loc[item, 'mafF_a'] = magF_a
                df.loc[item, 'BEXT'] = magF_a/magF_u
            except:
                df.loc[item, 'mafF_a'] = ''
                df.loc[item, 'BEXT'] = -5

            # placeholder for future implementation of site difference check
            mag_sites_difference(item)


    df.to_csv(datalist.replace(".csv", '_out' + '.csv'))
    # sieve(datalist.replace('.csv', '_out.csv'), 'magF_u', min_mag_field)


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

# vasp_results_dir = 'D:/MCES/MP/batch1/outdir'

wdatadir_structure = 'D:/MCES/MP/datadir/'
wdatalist = 'D:/MCES/MP/step4_success_sieved.csv'

# datalist = 'X:/MCES/Aflow/datalist_updated_sieved.mag.field_sieved.mag.sites_no.duplicates_beforeRun_afterRun_success_sieved_out.csv'
outdir = 'D:/MCES/MP/outdir'
# outdir = ('MnAs_o')
# wdatadir_structure = 'datadir/'
vasp_results_dir = outdir
vasp_results_BEXT = 'D:/MCES/MP/outdir_BEXT_1/out'
# screener_before('X:/MCES/MP/step2.csv', 'MP')

# sieve_temp02(wdatalist)

# print(calculate_mag_field(16.1092, 107.51))

# sieve_temp01('X:/MCES/MP/initial attempt - apparently wrong data from MRESTER or they had an update/outdir/datalist_lattfix_updated_sieved.mag.field_sieved.mag.sites_no.duplicates_beforeRun_afterRun_success_sieved_complete.csv')

# screener_after('X:/MCES/MP/initial attempt - apparently wrong data from MRESTER or they had an update/outdir/datalist_lattfix_updated_sieved.mag.field_sieved.mag.sites_no.duplicates_beforeRun_afterRun_success_sieved_complete_no_problem.csv')


# screener_after(wdatalist)

# sieve('D:/MCES/datalist_updated_sieved.mag.field_sieved.mag.sites_no.duplicates_beforeRun_afterRun_success_sieved_out.csv', 'magF_u', 0.45)

# magnetoelastic_max('mp-1200', 'a', plot=True)

# sieve_temp01('D:/MCES//MP/step3_failed_sieved_accuracy.csv')
# do_datalist('D:/MCES//MP/step3_failed_sieved_accuracy_dE.csv')
# duplicates('C:/Users/vikva/YandexDisk/Work/PycharmProjects/MCEScreening/Displayer/step3_success_sieved_out_and_extra.csv')
# rare_earth_content('mp-1208661')
screener_after('step4_success_sieved.csv')
# do_datalist('C:/Users/vikva/YandexDisk/Work/PycharmProjects/MCEScreening/Displayer/step4_success_sieved_out.csv')

# print(calculate_mag_field(9.288, 184.50))