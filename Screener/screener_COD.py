import pandas as pd
import numpy as np
import tqdm
from pymatgen.io.cif import CifParser
from pymatgen.symmetry.analyzer import SpacegroupAnalyzer
from pymatgen.core import Composition

### LOCAL IMPORTS ###
from Tools import POSCAR_reader
#from Screener import bezier
import time

from collections import Counter
import threading
import _thread   # import thread in python2

from gemmi import cif
import re

class timeout():
    def __init__(self, time):
        self.time= time
        self.exit=False

    def __enter__(self):
        threading.Thread(target=self.callme).start()

    def callme(self):
        time.sleep(self.time)
        if self.exit==False:
            _thread.interrupt_main()  # use thread instead of _thread in python2

    def __exit__(self, a, b, c):
        self.exit=True


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


def composition_analysis(ID):
    """ """
    rare_earths = ['Sc', 'Y', 'La', 'Ce', 'Pr', 'Nd', 'Pm', 'Sm', 'Eu', 'Gd', 'Tb', 'Dy', 'Ho', 'Er', 'Tm', 'Yb', 'Lu']
    magnetic = ['Ti', 'V', 'Cr', 'Mn', 'Fe', 'Co', 'Ni',        # Cu, Sc, Y, and everything else is considered not magnetic
                 'Nb', 'Mo', 'Ru', 'Rh', 'Pd',                  # we sieve out Cd at earlier step, so I could have omitted it from this list but this approach is more "universal"
                 'Ce', 'Pr', 'Nd', 'Sm', 'Eu', 'Gd', 'Tb', 'Dy', 'Ho', 'Er', 'Tm', 'Yb']
    oxygen = ['O']
    rare_total = 0
    oxygen_total = 0
    active_total = 0

    eltable = np.empty([0, 2])

    structure_file = wdatadir_structure + str(ID) + '.cif'
    structure = (CifParser(structure_file, occupancy_tolerance=0.9).get_structures())[0]

    a = SpacegroupAnalyzer(structure, symprec=1e-4)
    lat_type = a.get_crystal_system()
    if lat_type == 'trigonal':
        lat_type = 'rhombohedral'

    space_group_pymatgen = a.get_space_group_number()

    number = len(structure.cart_coords)
    el_list0 = [str(i).replace('Element ', '') for i in structure.species]
    el_list = [re.sub(r'[^a-zA-Z\s]', u'', str(j), flags=re.UNICODE) for j in el_list0]

    c = Counter(el_list)
    for i in c.keys():
        if i in oxygen:
            oxygen_total = oxygen_total + c[i]
        if i in rare_earths:
            rare_total = rare_total + c[i]
        if i in magnetic:
            active_total = active_total + c[i]
    o_content = 100*oxygen_total/number
    rare_content = 100*rare_total/number
    mag_active_content = 100*active_total/number

    comp = Composition(structure.composition)
    price_estimate = 0.0
    el_count = len(comp.elements)

    for i in comp.elements:
        element = re.sub(r'[^a-zA-Z\s]', u'', str(i), flags=re.UNICODE)
        price_per_kg = price_list.loc[str(element), 'price EUR/kg']
        if price_per_kg == '-':
            price_estimate = 'unknown'
            break
        else:
            price_el = comp.get_wt_fraction(i) * float(price_per_kg)
            price_estimate = price_estimate + price_el

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
    return number, rare_content, o_content, mag_active_content, price_estimate, site_counter, lat_type, space_group_pymatgen, el_count


def problematic_cif_parser_composition(datalist):
    counter = 0
    df = pd.read_csv(datalist, index_col=0, sep=',')
    with tqdm.tqdm(total=len(df.index)) as pbar:        # A wrapper that creates nice progress bar
        pbar.set_description("Processing datalist")
        for item in df.index.tolist():
            pbar.update(1)                              # Updating progress bar at each step

            try:
                structure_file_path = wdatadir_structure + str(item) + '.cif'
                structure_file = cif.read_file(structure_file_path)
                block = structure_file.sole_block()

                chemical_formula = block.find_value('_chemical_formula_sum')
                tot_number = len(block.find_loop('_atom_site_label'))

                eltable = np.empty(0, dtype=([('el', '<U5'), ('num', '<f8')]))
                elements_numbers = chemical_formula.strip("'").split(' ')
                for i in elements_numbers:
                    match = re.compile(r'[A-Za-z]+|-?\d+\.\d+|\d+|\W')
                    pair = match.findall(i)
                    el = str(pair[0])
                    if len(pair) > 1:
                        if pair[1] == '.':
                            num = float(pair[2])/10
                        else:
                            num = (float(pair[1]))
                    else:
                        num = 1
                    eltable = np.append(eltable, np.array([(el, num)], dtype=eltable.dtype), axis=0)

                number = eltable['num'].sum()

                magnetic = ['Ti', 'V', 'Cr', 'Mn', 'Fe', 'Co', 'Ni',  # Cu, Sc, Y, and everything else is considered not magnetic
                            'Nb', 'Mo', 'Ru', 'Rh', 'Pd',  # we sieve out Cd at earlier step, so I could have omitted it from this list but this approach is more "universal"
                            'Ce', 'Pr', 'Nd', 'Sm', 'Eu', 'Gd', 'Tb', 'Dy', 'Ho', 'Er', 'Tm', 'Yb']
                np.in1d(eltable['el'], magnetic)
                active_total = (eltable[np.in1d(eltable["el"], magnetic)]["num"].sum())
                mag_active_content = 100 * active_total / number

                df.loc[item, 'formula'] = chemical_formula
                df.loc[item, 'number_elements'] = number
                df.loc[item, 'number_of_atoms'] = int(tot_number)
                df.at[item, 'mag_active'] = round(mag_active_content, 2)
            except:
                counter = counter + 1
                print('cannot read cif file for ', item)
    df.to_csv(datalist.replace('.csv', '_updated.csv'))


def problematic_cif_parser_magsites(datalist):
    counter = 0
    df = pd.read_csv(datalist, index_col=0, sep=',')
    with tqdm.tqdm(total=len(df.index)) as pbar:  # A wrapper that creates nice progress bar
        pbar.set_description("Processing datalist")
        for item in df.index.tolist():
            pbar.update(1)  # Updating progress bar at each step

            try:
                magnetic = ['Ti', 'V', 'Cr', 'Mn', 'Fe', 'Co', 'Ni',  # Cu, Sc, Y, and everything else is considered not magnetic
                            'Nb', 'Mo', 'Ru', 'Rh', 'Pd',  # we sieve out Cd at earlier step, so I could have omitted it from this list but this approach is more "universal"
                            'Ce', 'Pr', 'Nd', 'Sm', 'Eu', 'Gd', 'Tb', 'Dy', 'Ho', 'Er', 'Tm', 'Yb']

                eltable = np.empty([0, 2])

                structure_file = wdatadir_structure + str(item) + '.cif'
                structure = (CifParser(structure_file, occupancy_tolerance=0.9).get_structures())[0]

                analyzed_structure = str(SpacegroupAnalyzer(structure).get_symmetrized_structure()).split('\n')
                for num, line in enumerate(analyzed_structure):
                    if num > 7:
                        elem = re.sub(r'[^a-zA-Z\s]', u'', str(line.split()[1]), flags=re.UNICODE)  # element symbol is after x y z coordinates, we also get rid of everything that is not a letter
                        wyckoff = line.split()[-1]  # wyckoff position is the last word in the string
                        eltable = np.append(eltable, [[elem, wyckoff]], axis=0)
                mask = np.in1d(eltable[:, 0], magnetic)  # Check what entries in Wyckoff list are in list of magnetic atoms
                eltable2 = eltable[mask]  # new array of only magnetic ones using previous mask
                unique_keys, mask2 = np.unique(eltable2[:, 1],
                                               return_index=True)  # check what entries of the new array have unique wyckoff sites - creates a list of indicies corresponding to said entries
                mag_sites = np.size(eltable2[mask2], axis=0)  # get number of unique sites counting number of obtained indicies
                df.loc[item, 'mag_sites'] = mag_sites  # mag_sites_calculator_MP(item)

            except:
                counter = counter + 1
                # print('cannot read cif file for ', item)
    df.to_csv(datalist.replace('.csv', '_updated.csv'))


def screener_before(datalist, database_type):
    # I feel that performance may not be optimal - making two loops does not seem reasonable
    # but I have to test if working with two df simultaneously is faster...
    """First main function that works with screening database. Used to apply criteria that don't require calculations.
    Loops over all database entries and applies initial screening parameters/checks.
    Creates a shorter database that we will use for submitting calculations"""

    # Criteria limits
    min_site_number = 1
    min_mag_field = 0.45
    min_active = 50  # in %
    counter = 0
    df = pd.read_csv(datalist, index_col=0, sep=',')
    fail_file = open(fail_file_path, "w")
    fail_file.write('ID,\n')
    with tqdm.tqdm(total=len(df.index)) as pbar:        # A wrapper that creates nice progress bar
        pbar.set_description("Processing datalist")
        for item in df.index.tolist():
            pbar.update(1)                              # Updating progress bar at each step

            try:
                with timeout(30):
                    number_atoms, rare_content, oxygen_content, mag_active_content, price, mag_sites, lat_type, space_group_pymatgen, el_count = composition_analysis(item)
                    df.at[item, 'number_atoms'] = number_atoms
                    # Rare earth content into datalist
                    df.at[item, 'rare_content'] = round(rare_content, 0)
                    # Oxygen content into datalist
                    df.at[item, 'oxygen_content'] = round(oxygen_content, 0)
                    # Write percentage of magnetically active elements
                    df.at[item, 'mag_active'] = round(mag_active_content, 2)
                    # Price estimate
                    df.at[item, 'price'] = round(price, 2)
                    # number of sites into datalist:
                    df.at[item, 'mag_sites'] = mag_sites # mag_sites_calculator_MP(item)
                    # symmetry type
                    df.loc[item, 'lattice_system'] = lat_type
                    # space group detected by pymatgen
                    df.loc[item, 'space_group_pymatgen'] = space_group_pymatgen
                    # number of distinct elements
                    df.loc[item, 'el_count'] = el_count

            except:
                counter = counter + 1
                fail_file.write(str(item) + ',\n')
    fail_file.close()
    # Wrire an updated datalist back to file
    print(counter)
    # df2.to_csv(datalist.replace(".csv", '_failed_cif_out' + '.csv'))
    df.to_csv(datalist.replace('.csv', '_updated.csv'))

    # drops everything that does not fit the criteria, creating new datalist file at each sieve
    duplicates(datalist.replace('.csv', '_updated.csv'))
    # sieve(datalist.replace('.csv', '_updated.csv'), 'mag_sites', min_site_number)
    sieve(datalist.replace('.csv', '_updated_no.duplicates.csv'), 'mag_active', min_active)
    if database_type == 'MP' or database_type == 'aflow':
        sieve(datalist.replace('.csv', '_updated_sieved.sieved.mag.sites_sieved.mag.active.csv'), 'mag_field', min_mag_field)

    print("Done")


def sieve(datalist, sieve_type, sieve_size):
    """Removes entries from datalist according to selected sieving criteria"""

    print("\nSieving by " + str(sieve_type) + " with cutoff set as " + str(sieve_size))
    df = pd.read_csv(datalist, index_col=0, sep=',')
    # for item in df.index.tolist():
    #     if df.loc[item, sieve_type] <= sieve_size:
    #         df = df.drop([item], axis=0)
    # df[[sieve_type]] = df[[sieve_type]].apply(pd.to_numeric)
    df = df[df[sieve_type] > sieve_size]
    df.to_csv(datalist.replace(".csv", '_sieved.' + sieve_type.replace('_', '.') + '.csv'))


# ------------------------------------------------------------------------------------------------------- #
#  _____                                           _       _____ _             _     _   _                #
# /  __ \                                         | |     /  ___| |           | |   | | | |               #
# | /  \/ ___  _ __ ___  _ __ ___   __ _ _ __   __| |___  \ `--.| |_ __ _ _ __| |_  | |_| | ___ _ __ ___  #
# | |    / _ \| '_ ` _ \| '_ ` _ \ / _` | '_ \ / _` / __|  `--. \ __/ _` | '__| __| |  _  |/ _ \ '__/ _ \ #
# | \__/\ (_) | | | | | | | | | | | (_| | | | | (_| \__ \ /\__/ / || (_| | |  | |_  | | | |  __/ | |  __/ #
#  \____/\___/|_| |_| |_|_| |_| |_|\__,_|_| |_|\__,_|___/ \____/ \__\__,_|_|   \__| \_| |_/\___|_|  \___| #
#                                                                                                         #
# ------------------------------------------------------------------------------------------------------- #
price_list = pd.read_csv('prices.csv', index_col=0, sep=',')  # read the prices file
### Setting which database we work with ###

# Aflow
# wdatadir_structure = '../Database/aflow/datadir_structure_relaxed/'
# wdatalist = '../Database/aflow/datalist_updated_sieved.mag.field_sieved.mag.sites.csv'

# Materials Project
# calculations = {'undeformed': '', 'Applied_Field': '', 'uniaxial': ['0.84', '0.9', '0.95', '1.05', '1.1', '1.16'], 'volumetric': ['0.95', '1.05']}
# wdatadir_structure = 'D:/MCES/MP/datadir/'
# wdatalist = 'D:/MCES/MP/outdir.csv'
# vasp_results_dir = 'D:/MCES/MP/outdir'


# COD
# fail_file_path = 'D:/MCES/COD/failed_cif_out_3rd.csv'
# wdatadir_structure = 'D:/MCES/COD/datadir/'
# wdatalist = 'D:/MCES/COD/datalist_COD.csv'


# ICSD
# fail_file_path = 'D:/MCES/ICSD/failed_cif_out_1st.csv'
# wdatadir_structure = 'D:/MCES/ICSD/datadir/'
# wdatalist = 'D:/MCES/ICSD/datalist.csv'
# screener_before(wdatalist, 'ICSD')

fail_file_path = 'D:/MCES/TESTS/Bocarsly/failed_cif_out_1st.csv'
wdatalist = 'D:/MCES/TESTS/Bocarsly/datalist.csv'
wdatadir = 'D:/MCES/TESTS/Bocarsly/datadir/'
output_path = 'D:/MCES/TESTS/Bocarsly/inputdir/'
screener_before(wdatalist, 'COD')


# screener_after(wdatalist)
# sieve(wdatalist.replace('.csv', '_updated.csv'), 'mag_active', 33)
# screener_before(wdatalist, 'COD')
# # screener_before('D:/MCES/COD/failed_cif_out_sieved.mag.active_sieved.mag.sites.csv', 'COD')
# # problematic_cif_parser_magsites('D:/MCES/COD/failed_cif_out_sieved.mag.active.csv')
# # cif_parser(fail_file_path)
# # cif_parser('D:/MCES/COD/failed_cif_out1.csv')
