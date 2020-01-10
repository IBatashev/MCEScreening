def read(structure_file):
    """Takes structure file as input and returns section containing POSCAR as list of strings"""

    poscar_content = []
    with open(structure_file, 'r') as f:  # need to get POSCAR content from big structure file
        for line in f:
            if 'SPRIM' in line:  # subsection we are interested in goes after line with word "Representative"
                for line in f:  # now you are at the lines you want
                    if 'SCONV' in line:  # Ends before line with "WYCCAR"
                        break
                    else:
                        poscar_content.append(line)
    return poscar_content

 #
 #
 # with open(structure_file, 'r') as f:    # need to get all atoms from Wyckoff subsection of structure file
 #        for line in f:
 #            if 'Representative' in line:    # subsection we are interested in goes after line with word "Representative"
 #                for line in f:              # now you are at the lines you want
 #                    if 'WYCCAR' in line:    # Ends before line with "WYCCAR"
 #                        break
 #                    else:
 #                        elem = line.split()[3]  # element symbol is after Wyckoff x y z coordinates so 4th item in list
 #                        if elem == 'A':         # check if notation is shitty - in this case first element is represented with letter A (fortunately there is no element in periodic table labeled with A)
 #                            bad_notation = True
 #                        if bad_notation == True:# Apply correction for bad notation replacing all alphabet letters with proper chemical symbols taken from list of elements in composition
 #                            n = rename_list.index(elem)
 #                            temp_line = (df.loc[ID, 'species']).strip('[').strip(']').strip(' ')
 #                            temp_line2 = temp_line.replace("'", "")
 #                            species_list = np.asarray(temp_line2.split('; '))
 #                            elem = species_list[n]
 #                        elemlist = np.append(elemlist, elem)
 #    for k in elemlist:
 #        if k in magnetic:
 #            site_counter = site_counter + 1
 #    return site_counter