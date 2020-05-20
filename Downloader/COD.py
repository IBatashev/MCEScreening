import os
import tqdm
import numpy as np
import pandas as pd
from pymatgen.io.cif import CifParser
import pymatgen.symmetry.analyzer
from pymatgen.io.vasp import Poscar

# _cod_database_code
# _chemical_formula_sum
# _chemical_formula_structural
# _symmetry_cell_setting lattice system
# _space_group_IT_number spacegroup
# elements -split from  _chemical_formula_structural
# _cell_volume
# _journal_paper_doi
# _publ_section_title

# unique sites -calculate from structure and elements


def make_path_list():
    directory = 'D:/cif/'
    listing = [1, 2, 3, 4, 5, 6, 7, 8, 9]
    with open('all_paths', "w+") as f:
        for i in listing:
            listing2 = (os.listdir(directory+str(i)))
            for k in listing2:
                listing3 = (os.listdir(directory+str(i) + '/' + str(k)))
                for j in listing3:
                    listing4 = os.listdir(directory + str(i) + '/' + str(k) + '/' + str(j))
                    for l in listing4:
                        f.write(directory + str(i) + '/' + str(k) + '/' + str(j) + '/' + str(l))
                        f.write('\n')


def scan(pathlist):
    with open('datalist_COD.csv', 'a') as datalist:  # we create a csv file to write info into
        datalist.write("COD_ID,path,pretty_formula,compound,lattice_system,spacegroup,species,volume_cell,mag_sites,comment1,doi\n")
    no_pretty_formula_counter = 0
    total_compounds = len(open(pathlist).readlines())
    with tqdm.tqdm(total=total_compounds) as pbar:        # A wrapper that creates nice progress bar
        pbar.set_description("Processing")
        with open(pathlist, 'r') as f:  # need to get POSCAR content from big structure file
            for path in f:
                pbar.update(1)
                with open(path.strip('\n'), 'r') as cif:
                    content = cif.read()
                    try:
                        pretty_formula = content.split('_chemical_formula_sum')[1].split('\n')[0].strip(' ').strip("'").strip()
                    except:
                        pretty_formula = 'na'
                        no_pretty_formula_counter = no_pretty_formula_counter + 1
                        pass
                    try:
                        elemets = pretty_formula.translate({ord(i): None for i in '1234567890'}).split()
                    except:
                        elemets = 'na'
                        pass

                    necessary = ["Mn", "Fe", "Co", "Ni", "Cu"] # at least one of these must be present
                    banlist = ["Re", "Os", "Ir", "Pt", "Au", "In", "Tc",  # Expensive or Limited in supply
                               "Be", "As", "Cd", "Ba", "Hg", "Tl", "Pb", "Ac",  # Health Hazard
                               "Cs", "Pa", "Np", "U", "Pu", "Th",  # Radioactive
                               "He", "Ne", "Ar", "Kr", "Xe"]  # Noble gases

                    necessary_match = [i for i in necessary if i in elemets]
                    if necessary_match:
                        banlist_match = [i for i in banlist if i in elemets]

                        if not banlist_match:
                            try:
                                COD_id = content.split('_cod_database_code')[1].split('\n')[0].strip(' ')
                            except:
                                COD_id = 'na'
                            try:
                                compound = content.split('_chemical_formula_structural')[1].split('\n')[0].strip(' ').strip("'")
                            except:
                                compound = 'na'
                            try:
                                lattice_system = content.split('_symmetry_cell_setting')[1].split('\n')[0].strip(' ').strip("'")
                            except:
                                lattice_system = 'na'
                            try:
                                spacegroup = content.split('_space_group_IT_number')[1].split('\n')[0].strip(' ').strip("'")
                            except:
                                spacegroup = 'na'

                            try:
                                volume = content.split('_cell_volume')[1].split('\n')[0].strip(' ').strip("'")
                            except:
                                volume = 'na'
                            try:
                                doi = content.split('_journal_paper_doi')[1].split('\n')[0].strip(' ').strip("'")
                            except:
                                doi = 'na'
                            try:
                                comment1 = content.split('_publ_section_title')[1].split('\n')[2].split('_journal_name_full')[0]
                            except:
                                comment1 = 'na'

                            # Now we write all parsed results into .csv file

                            newrow = str(
                                COD_id) + ',' + str(
                                path.strip('D:/').strip('\n')) + ',' + str(
                                pretty_formula) + ',' + str(
                                compound) + ',' + str(
                                lattice_system) + ',' + str(
                                spacegroup) + ',' + str(
                                elemets).replace(',', ';') + ',' + str(
                                volume) + ',' + str(
                                0) + ',' + str(
                                comment1).replace(',', ';') + ',' + str(
                                doi).strip(',') + ',' + '\n'
                            with open('datalist_COD.csv', 'a') as datalist:
                                datalist.write(newrow)

    print('No chemical composition given in cif for ', no_pretty_formula_counter, ' entries')


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

    structure_file = wdatadir_structure + df.loc[str(ID), 'path']

    parser = CifParser(structure_file)  # create structure from .cif using pymatgen
    structure = parser.get_structures()[0]
    print(pymatgen.symmetry.analyzer.SpacegroupAnalyzer(structure))


    # p_content = Poscar(structure)  # create poscar from structure
    # p = str(p_content).split('\n')  # These two lines are to remodel poscar_content to the same
    # poscar_content = [string + '\n' for string in p]  # data format we get from aflow - a list of strings



    # with open(structure_file, 'r') as f:    # need to get all atoms from Wyckoff subsection of structure file
    #     for line in f:
    #         if 'Representative' in line:    # subsection we are interested in goes after line with word "Representative"
    #             for line in f:              # now you are at the lines you want
    #                 if 'WYCCAR' in line:    # Ends before line with "WYCCAR"
    #                     break
    #                 else:
    #                     elem = line.split()[3]  # element symbol is after Wyckoff x y z coordinates so 4th item in list
    #                     wyckoff = line.split()[4]+line.split()[5]
    #                     if elem == 'A':         # check if notation is shitty - in this case first element is represented with letter A (fortunately there is no element in periodic table labeled with A)
    #                         bad_notation = True
    #                     if bad_notation == True:# Apply correction for bad notation replacing all alphabet letters with proper chemical symbols taken from list of elements in composition
    #                         n = rename_list.index(elem)
    #                         temp_line = (df.loc[ID, 'species']).strip('[').strip(']').strip(' ')
    #                         temp_line2 = temp_line.replace("'", "")
    #                         species_list = np.asarray(temp_line2.split('; '))
    #                         elem = species_list[n]
    #                     eltable = np.append(eltable, [[elem, wyckoff]], axis=0)
    # # This section with masks probably could be done better and shorter but with such small arrays it shouldn't matter much

    # mask = np.in1d(eltable[:, 0], magnetic)  # Check what entries in Wyckoff list are in list of magnetic atoms
    # eltable2 = eltable[mask]                 # new array of only magnetic ones using previous mask
    # unique_keys, mask2 = np.unique(eltable2[:, 1], return_index=True)  # check what entries of the new array have unique wyckoff sites - creates a list of indicies corresponding to said entries
    # site_counter = np.size(eltable2[mask2], axis=0)                    # get number of unique sites counting number of obtained indicies
    # return site_counter


wdatalist = 'datalist_COD.csv'
wdatadir_structure = 'D:/'

df = pd.read_csv(wdatalist, index_col=0, sep=',')
# with tqdm.tqdm(total=len(df.index)) as pbar:  # A wrapper that creates nice progress bar
#     pbar.set_description("Processing datalist")
#     for item in df.index.tolist():
#         pbar.update(1)  # Updating progress bar at each step

mag_sites_calculator(1004018)
