import phonopy

def mag_sites_calculator(ID):
    """Determines how many unique magnetic sites are present in structure.
    Takes ID as imput, looks up poscar file in the datadir and returns number of sites as integer
    """
    poscar_file = '../Database/datadir/'+str(ID)
    ### List of Magnetic Atoms
    magnetic = ['Cr', 'Mn', 'Fe', 'Co', 'Ni', 'Gd'] # need to add all of them later

    ### Open POSCAR with phonopy and get Wyckoff letters and chemical symbols from it
    # maybe should try with ase????????
    ph = phonopy.load(unitcell_filename=poscar_file)
    sym_list = ph.get_symmetry().get_Wyckoff_letters()
    at_type_len = len(ph.unitcell.get_positions())
    at_type_list = []

    ### Open POSCAR and get element symbols perhaps make this a separrate function for both runner and screener...
    poscar = open(poscar_file, "r")
    list = []
    for line in poscar:
        list.append(line)
    poscar.close()

    for i in range(7, 7 + at_type_len): # POSCAR files from aflowlib have 7 lines of text before atomic coordinates! so we start at line 8
        l = str.split(list[i], '  ')    # Aflow POSCAR has symbols for elements listed after their coordinate
        at_type_list.append(str(l[4]))  # create a list of all elements in the structure

    ### Check how many unique MAGNETIC sites are present
    at_list = []
    for num, val in enumerate(sym_list):      ### build list of magnetic positions
        if at_type_list[num] in magnetic:
           at_list.append(val)
    num_unique = len(set(at_list))            ### check how many are unique
    return num_unique                         ### return resulting number of sites

    ### Another attempt to do same as above need to test both further

    # print(ph.get_symmetry().get_dataset())
    # for i, vall in enumerate(at_type_list):
    #     sym_point = ph.get_symmetry().get_site_point_group()
    #     # sym_list = sym_list.append(ph.get_symmetry().get_pointgroup(sym_point))
    #     print(sym_point)
    # print(sym_list)
    # print(at_type_list)