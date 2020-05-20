from pymatgen import MPRester
m = MPRester('LXELGKLFgVPOfRTSj')

out_path = 'MP_structures/'
filename = 'temp_file_cif'


def chose_ids_to_download():
    """Creates file(s) containing Materials Project ids matching first criteria that we can next supply to download()"""

    #  Does not contain any of these                                        because:
    banlist = ["Re", "Os", "Ir", "Pt", "Au", "In", "Tc",                    # Expensive or Limited in supply
               "Be", "As", "Cd", "Ba", "Hg", "Tl", "Pb", "Ac",              # Health Hazard
               "Cs", "Pa", "Np", "U", "Pu", "Th",                           # Radioactive
               "He", "Ne", "Ar", "Kr", "Xe"]                                # Noble gases

    data = m.query(criteria={"elements": {"$in": ["Mn", "Fe", "Co", "Ni", "Cu"],
                                          "$nin": banlist},
                             "magnetism.total_magnetization_normalized_vol": {"$gte": 0.03861276264262763},  # the number corresponds to internal field of 0.45T
                             "magnetism.num_unique_magnetic_sites": {"$gte": 2},
                             "formation_energy_per_atom": {"$lt": 0}
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
                                "exp.tags",                             # essentially a comment for this MP entry, information on this material derived from experimental databases such as the ICSD
                                "icsd_ids",                             # (ICSD) ids for structures that have been deemed to be structurally similar to this material based on pymatgen's StructureMatcher algorithm.
                                "doi",                                  # some seem to point to corresponding MP webpage, which is not bad but less useful than actual paper...
                                "warnings",                             # warnings associated with the material
                                "e_above_hull",                         # The calculated energy above the convex hull from the phase diagram. An indication of how stable a material is. A stable material is on the hull and has an e_above_hull of 0. A larger positive number indicates increasing instability.
                                "cif"]                                  # the structure in the CIF format.
                    )
    print('Downloaded information on ', len(rdata), 'compounds, now saving...')

    with open('datalist_MP.csv', 'a') as f:  # we create a csv file to write info into
        f.write("ID,pretty_formula,compound,energy_cell,energy_atom,lattice_system,spacegroup,species,volume_cell,moment_cell,mag_field,mag_sites_MP,mag_sites,mag_type,comment1,icsd_ids,doi,comment2,e_above_hull,magnetization_norm_vol\n")

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
            d["exp.tags"]).replace(',', ';') + ',' + str(
            d["icsd_ids"]).replace(',', ';') + ',' + str(
            d["doi"]) + ',' + str(
            d["warnings"]).replace(',', ';') + ',' + str(
            d["e_above_hull"]) + ',' + str(
            d["magnetism.total_magnetization_normalized_vol"]) + '\n'

        with open('datalist_MP.csv', 'a', encoding="utf-8") as f:
            f.write(newrow)

        with open(out_path + str(d['material_id']) + str('.cif'), "w+", encoding="utf-8") as f:
            f.write(d['cif'])
    print('Done')


# chose_ids_to_download()
download('ids_to_download_all')


# https://github.com/materialsproject/mapidoc

# CITING
# Ong, S. P.; Cholia, S.; Jain, A.; Brafman, M.; Gunter, D.; Ceder, G.;
# Persson, K. a. The Materials Application Programming Interface (API): A
# simple, flexible and efficient API for materials data based on
# REpresentational State Transfer (REST) principles, Comput. Mater. Sci.,
# 2015, 97, 209–215. doi:10.1016/j.commatsci.2014.10.037.
