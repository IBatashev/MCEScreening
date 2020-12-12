import pandas as pd
import tqdm
from pymatgen.io.vasp import Poscar
from Tools import POSCAR_reader
from pymatgen import MPRester
m = MPRester('LXELGKLFgVPOfRTSj')

# https://github.com/materialsproject/mapidoc

# CITING
# Ong, S. P.; Cholia, S.; Jain, A.; Brafman, M.; Gunter, D.; Ceder, G.;
# Persson, K. a. The Materials Application Programming Interface (API): A
# simple, flexible and efficient API for materials data based on
# REpresentational State Transfer (REST) principles, Comput. Mater. Sci.,
# 2015, 97, 209–215. doi:10.1016/j.commatsci.2014.10.037.


def chose_ids_to_download():
    """Creates file(s) containing Materials Project ids matching first criteria that we can next supply to download()"""

    #  Does not contain any of these                                        because:
    banlist = ["Re", "Os", "Ir", "Pt", "Au", "In", "Tc",                    # Expensive or Limited in supply
               "Be", "As", "Cd", "Ba", "Hg", "Tl", "Pb", "Ac",              # Health Hazard
               "Cs", "Pa", "Np", "U", "Pu", "Th",                           # Radioactive
               "He", "Ne", "Ar", "Kr", "Xe"]                                # Noble gases

    data = m.query(criteria={"elements": {"$in": ["Mn", "Fe", "Co", "Ni", "Cr"],
                                          "$nin": banlist},
                             "total_magnetization": {"$gt": 0},
                             "formation_energy_per_atom": {"$lt": 0},
                             "magnetism.total_magnetization_normalized_vol": {"$gte": 0.03861276264262763},  # the number corresponds to internal field of 0.45T
                             "magnetism.num_unique_magnetic_sites": {"$gte": 2}
                             },
                   properties=["material_id"])

    print('Found ', len(data), 'matching compounds.')

    with open('ids_to_download_all', "w+") as f:
        for num, d in enumerate(data):
            f.write((d["material_id"]))
            f.write("\n")

    # batch splitter - to create sublists containing MP_id to download so that we can safely do actual downloading in several batches
    batch_filename_list = []
    lines_per_file = 2000
    smallfile = None
    count = 0
    with open('ids_to_download_all') as bigfile:
        for lineno, line in enumerate(bigfile):
            if lineno % lines_per_file == 0:
                if smallfile:
                    smallfile.close()
                count = count + 1
                small_filename = 'ids_to_download_part_{}'.format(count)
                batch_filename_list.append(small_filename)
                smallfile = open(small_filename, "w")
            smallfile.write(line)
        if smallfile:
            smallfile.close()
    print('Created ', count, ' sublists to download:\n', batch_filename_list)
    with open('ids_to_download_list_of_sublists', "w+") as listfile:
        listfile.write(str(batch_filename_list).strip('[').strip(']'))


def download(file_with_ids_to_download):
    """Takes a file containing a column of Materials Project ids and downloads all
     relevant info into a .csv file and a bunch of cif files"""

    downloadlist = []
    with open(file_with_ids_to_download, 'r') as f:
        for line in f:
            downloadlist.append(line.strip("\n"))
    # LIST of what to save for more details/options see https://github.com/materialsproject/mapidoc/tree/master/materials
    # Note that the data returned is always a list of dicts.
    rdata = m.query(criteria={"material_id": {"$in": downloadlist}},
                    properties=["material_id",                          # essentially a MP Id BUT each compound in MP has multiple calculations each with such number...
                                "pretty_formula",                       # formula where the element amounts are normalized. E.g., "Li2O"
                                "full_formula",                         # full explicit formula for the unit cell, e.g., "Lu2Al4"
                                "final_energy_per_atom",                #
                                "final_energy",                         # Calculated vasp energy for structure
                                "spacegroup.crystal_system",            #
                                "spacegroup.number",                    #
                                "elements",                             # separate elements
                                "volume",                               # The volume of the unit cell of the material used in the calculation. Need tocheck if it corresponds to cif volume, if not, cif is more reliable
                                "magnetism.total_magnetization",        #
                                "magnetism.num_unique_magnetic_sites",
                                "magnetism.total_magnetization_normalized_vol",
                                "magnetism.ordering",                   # type of magnetic order
                                "magnetism.magmoms",
                                "exp.tags",                             # essentially a comment for this MP entry, information on this material derived from experimental databases such as the ICSD
                                "icsd_ids",                             # (ICSD) ids for structures that have been deemed to be structurally similar to this material based on pymatgen's StructureMatcher algorithm.
                                "doi",                                  # some seem to point to corresponding MP webpage, which is not bad but less useful than actual paper...
                                "warnings",                             # warnings associated with the material
                                "e_above_hull",                         # The calculated energy above the convex hull from the phase diagram. An indication of how stable a material is. A stable material is on the hull and has an e_above_hull of 0. A larger positive number indicates increasing instability.
                                # "cif"]                                # the structure in the CIF format.
                                "cifs.conventional_standard"]
                    )
    print('Downloaded information on ', len(rdata), 'compounds, now saving...')

    with open('datalist.csv', 'a') as f:  # we create a csv file to write info into
        f.write("ID,pretty_formula,compound,energy_cell,energy_atom,lattice_system,spacegroup,species,volume_cell,moment_cell,mag_field,mag_sites_MP,mag_sites,mag_type,magmom,comment1,icsd_ids,doi,comment2,e_above_hull,magnetization_norm_vol\n")

    for num, d in enumerate(rdata):
        newrow = str(
            d["material_id"]) + ',' + str(
            d["pretty_formula"]) + ',' + str(
            d["full_formula"]) + ',' + str(
            d["final_energy"]) + ',' + str(
            d["final_energy_per_atom"]) + ',' + str(
            d["spacegroup.crystal_system"]) + ',' + str(
            d["spacegroup.number"]) + ',' + str(
            d["elements"]).replace(',', ';') + ',' + str(
            d["volume"]) + ',' + str(
            d["magnetism.total_magnetization"]) + ',' + str(
            0) + ',' + str(
            d["magnetism.num_unique_magnetic_sites"]) + ',' + str(
            0) + ',' + str(
            d["magnetism.ordering"]) + ',' + str(
            d["magnetism.magmoms"]).replace(',', ';') + ',' + str(
            d["exp.tags"]).replace(',', ';') + ',' + str(
            d["icsd_ids"]).replace(',', ';') + ',' + str(
            d["doi"]) + ',' + str(
            d["warnings"]).replace(',', ';') + ',' + str(
            d["e_above_hull"]) + ',' + str(
            d["magnetism.total_magnetization_normalized_vol"]) + '\n'

        with open('datalist.csv', 'a', encoding="utf-8") as f:
            f.write(newrow)

        with open('datadir/' + str(d['material_id']) + str('.cif'), "w+", encoding="utf-8") as f:
            # f.write(d['cif'])
            f.write(d['cifs.conventional_standard'])
    print('Done')


def lattice_type_fix(datalist, datadir):
    """MP has a trigonal lattice type which can be both HEX and RHL we want to distinguish them for our deformations
    This short script does this based on the space group"""

    prec = 0.1
    angle_prec = 5.0

    df = pd.read_csv(datalist, index_col=0, sep=',')
    with tqdm.tqdm(total=len(df.index)) as pbar:        # A wrapper that creates nice progress bar
        pbar.set_description("Processing datalist")

        for item in df.index.tolist():
            pbar.update(1)                              # Updating progress bar at each step
            structure_file = datadir + str(item) + '.cif'
            poscar_content = POSCAR_reader.read(structure_file)
            poscar_string = ''.join(poscar_content)  # merging poscar content in a single string so pymatgen can read it
            lat_type = df.loc[item, 'lattice_system']

            if lat_type == 'trigonal':
                poscar = Poscar.from_string(poscar_string)  # using pymatgen to acquire our structure from poscar content
                structure = poscar.structure
                group = ((structure.get_space_group_info(symprec=prec, angle_tolerance=angle_prec))[0])
                if 'R' in str(group):
                    df.loc[item, 'lattice_system'] = 'rhombohedral'
                else:
                    df.loc[item, 'lattice_system'] = 'hexagonal'
        df.to_csv(datalist.replace(".csv", '_lattfix' + '.csv'))


def id_download(id):
    """Displays chosen information for a single id"""

    rdata = m.query(criteria={"material_id": id},
                    properties=["material_id",                          # essentially a MP Id BUT each compound in MP has multiple calculations each with such number...
                                "task_id",
                                "pretty_formula",                       # formula where the element amounts are normalized. E.g., "Li2O"
                                "full_formula",                         # full explicit formula for the unit cell, e.g., "Lu2Al4"
                                "final_energy_per_atom",                #
                                "final_energy",                         # Calculated vasp energy for structure
                                "spacegroup.crystal_system",            #
                                "spacegroup.number",                    #
                                "elements",                             # separate elements
                                "volume",                               # The volume of the unit cell of the material used in the calculation. Need tocheck if it corresponds to cif volume, if not, cif is more reliable
                                "magnetism",        #
                                "magnetism.num_unique_magnetic_sites",
                                "magnetism.total_magnetization_normalized_vol",
                                "magnetism.ordering",                   # type of magnetic order
                                "magnetism.magmoms",
                                "exp.tags",                             # essentially a comment for this MP entry, information on this material derived from experimental databases such as the ICSD
                                "icsd_ids",                             # (ICSD) ids for structures that have been deemed to be structurally similar to this material based on pymatgen's StructureMatcher algorithm.
                                "doi",                                  # some seem to point to corresponding MP webpage, which is not bad but less useful than actual paper...
                                "warnings",                             # warnings associated with the material
                                "e_above_hull",                         # The calculated energy above the convex hull from the phase diagram. An indication of how stable a material is. A stable material is on the hull and has an e_above_hull of 0. A larger positive number indicates increasing instability.
                                # "cif"]                                # the structure in the CIF format.
                                "cifs.conventional_standard",
                                "input"]
                    )
    for num, d in enumerate(rdata):
        print(d["magnetism"])
        print(d["material_id"])
        # print(d["task_id"])

        print(d["volume"])


def download_custom(file_with_ids_to_download):
    """Takes a file containing a column of Materials Project ids and downloads chosen info into a .csv file
    Functionally an exact copy of download().
    This one is to be used for specific checks. To avoid changing the working one.
    For regular download one should use simple download()"""

    downloadlist = []
    with open(file_with_ids_to_download, 'r') as f:
        for line in f:
            downloadlist.append(line.strip("\n"))

    rdata = m.query(criteria={"material_id": {"$in": downloadlist}},
                    properties=["material_id",                          # essentially a MP Id BUT each compound in MP has multiple calculations each with such number...
                                "pretty_formula",                       # formula where the element amounts are normalized. E.g., "Li2O"
                                "full_formula",                         # full explicit formula for the unit cell, e.g., "Lu2Al4"
                                "run_type",
                                "is_hubbard",
                                "last_updated"]
                    )
    print('Downloaded information on ', len(rdata), 'compounds, now saving...')

    with open('datalist.csv', 'a') as f:  # we create a csv file to write info into
        f.write("ID,pretty_formula,compound,run_type,last_updated\n")

    for num, d in enumerate(rdata):
        newrow = str(
            d["material_id"]) + ',' + str(
            d["pretty_formula"]) + ',' + str(
            d["full_formula"]) + ',' + str(
            d["run_type"]) + ',' + str(
            d["is_hubbard"]) + ',' + str(
            d["last_updated"]) + '\n'

        with open('datalist.csv', 'a', encoding="utf-8") as f:
            f.write(newrow)

    print('Done')

# ------------------------------------------------------------------------------------------------------- #
#  _____                                           _       _____ _             _     _   _                #
# /  __ \                                         | |     /  ___| |           | |   | | | |               #
# | /  \/ ___  _ __ ___  _ __ ___   __ _ _ __   __| |___  \ `--.| |_ __ _ _ __| |_  | |_| | ___ _ __ ___  #
# | |    / _ \| '_ ` _ \| '_ ` _ \ / _` | '_ \ / _` / __|  `--. \ __/ _` | '__| __| |  _  |/ _ \ '__/ _ \ #
# | \__/\ (_) | | | | | | | | | | | (_| | | | | (_| \__ \ /\__/ / || (_| | |  | |_  | | | |  __/ | |  __/ #
#  \____/\___/|_| |_| |_|_| |_| |_|\__,_|_| |_|\__,_|___/ \____/ \__\__,_|_|   \__| \_| |_/\___|_|  \___| #
#                                                                                                         #
# ------------------------------------------------------------------------------------------------------- #

wdatadir = 'datadir/'
wdatalist = 'datalist.csv'

# # Routine:
# # step 1:
# chose_ids_to_download()
# # step 2:
# download('ids_to_download_all')
# # step 3:
# lattice_type_fix(wdatalist, wdatadir)

# Extra stuff:

# download_custom('ids')
download('ids')
# id_download('mp-778')