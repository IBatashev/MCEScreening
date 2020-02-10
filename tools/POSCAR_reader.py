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

