import pandas as pd
import os
import numpy as np
import tqdm
import re


def composition_analysis(formula):
    """ """
    rare_earths = ['Sc', 'Y', 'La', 'Ce', 'Pr', 'Nd', 'Pm', 'Sm', 'Eu', 'Gd', 'Tb', 'Dy', 'Ho', 'Er', 'Tm', 'Yb', 'Lu']
    magnetic = ['Ti', 'V', 'Cr', 'Mn', 'Fe', 'Co', 'Ni',        # Cu, Sc, Y, and everything else is considered not magnetic
                 'Nb', 'Mo', 'Ru', 'Rh', 'Pd',                  # we sieve out Cd at earlier step, so I could have omitted it from this list but this approach is more "universal"
                 'Ce', 'Pr', 'Nd', 'Sm', 'Eu', 'Gd', 'Tb', 'Dy', 'Ho', 'Er', 'Tm', 'Yb']
    oxygen = ['O']
    rare_total = 0
    oxygen_total = 0
    active_total = 0
    # eltable = np.empty([0, 2], dtype=[('el', '<U5'), ('num', '<f8')])

    eltable = np.empty(0, dtype=([('el', '<U5'), ('num', '<f8')]))
    # dtype = (('<U5, float64'))

    elements = re.sub(r'\d+', '', formula).replace('.', '').split(' ')
    elements = str(elements).replace(',', ';')

    elements_numbers = formula.split(' ')
    for i in elements_numbers:
        # match = re.match(r"([a-z]+)([0-9]+)", i, re.I)
        match = re.compile(r'[A-Za-z]+|-?\d+\.\d+|\d+|\W')
        pair = match.findall(i)
        el = str(pair[0])
        num = (float(pair[1]))
        eltable = np.append(eltable, np.array([(el, num)], dtype=eltable.dtype), axis=0)

    number = eltable['num'].sum()

    magnetic = ['Ti', 'V', 'Cr', 'Mn', 'Fe', 'Co', 'Ni',        # Cu, Sc, Y, and everything else is considered not magnetic
                 'Nb', 'Mo', 'Ru', 'Rh', 'Pd',                  # we sieve out Cd at earlier step, so I could have omitted it from this list but this approach is more "universal"
                 'Ce', 'Pr', 'Nd', 'Sm', 'Eu', 'Gd', 'Tb', 'Dy', 'Ho', 'Er', 'Tm', 'Yb']

    np.in1d(eltable['el'], magnetic)
    # # array([False, True, True, True, True], dtype=bool)

    active_total = (eltable[np.in1d(eltable["el"], magnetic)]["num"].sum())
    mag_active_content = 100 * active_total / number

    # el_list = [str(i).replace('Element ', '') for i in structure.species]
    #
    # c = Counter(el_list)
    # for i in c.keys():
    #     if i in oxygen:
    #         oxygen_total = oxygen_total + c[i]
    #     if i in rare_earths:
    #         rare_total = rare_total + c[i]
    #     if i in magnetic:
    #         active_total = active_total + c[i]
    # o_content = 100*oxygen_total/number
    # rare_content = 100*rare_total/number
    # mag_active_content = 100*active_total/number

    return elements, mag_active_content #rare_content, o_content, mag_active_content


def screener_before(datalist, database_type='ICSD'):
    """First main function that works with screening database. Used to apply criteria that don't require calculations.
    Loops over all database entries and applies initial screening parameters/checks.
    Creates a shorter database that we will use for submitting calculations"""

    # Criteria limits
    min_site_number = 1
    min_mag_field = 0.45
    min_active = 33  # in %
    counter = 0
    df = pd.read_csv(datalist, index_col=0, sep=',')
    with tqdm.tqdm(total=len(df.index)) as pbar:        # A wrapper that creates nice progress bar
        pbar.set_description("Processing datalist")
        for item in df.index.tolist():
            pbar.update(1)                              # Updating progress bar at each step
            #rare_content, oxygen_content, mag_active_content = composition_analysis(df.loc[item, 'formula'])
            el_list, mag_active_content = (composition_analysis(df.loc[item, 'pretty_formula']))
            df.at[item, 'species'] = el_list
            df.at[item, 'mag_active'] = round(mag_active_content, 2)

    df.to_csv(datalist.replace('.csv', '_updated.csv'))

# ------------------------------------------------------------------------------------------------------- #
#  _____                                           _       _____ _             _     _   _                #
# /  __ \                                         | |     /  ___| |           | |   | | | |               #
# | /  \/ ___  _ __ ___  _ __ ___   __ _ _ __   __| |___  \ `--.| |_ __ _ _ __| |_  | |_| | ___ _ __ ___  #
# | |    / _ \| '_ ` _ \| '_ ` _ \ / _` | '_ \ / _` / __|  `--. \ __/ _` | '__| __| |  _  |/ _ \ '__/ _ \ #
# | \__/\ (_) | | | | | | | | | | | (_| | | | | (_| \__ \ /\__/ / || (_| | |  | |_  | | | |  __/ | |  __/ #
#  \____/\___/|_| |_| |_|_| |_| |_|\__,_|_| |_|\__,_|___/ \____/ \__\__,_|_|   \__| \_| |_/\___|_|  \___| #
#                                                                                                         #
# ------------------------------------------------------------------------------------------------------- #

# wdatadir_structure = 'D:/MCES/COD/datadir/'
wdatalist = 'D:/MCES/ICSD/id_list2_no_O.csv'


screener_before(wdatalist)
