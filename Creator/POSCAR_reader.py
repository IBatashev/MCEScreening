from pymatgen.io.cif import CifParser
from pymatgen.io.vasp import Poscar


def read(structure_file):
    """Takes structure file as input and returns section containing POSCAR as list of strings ending with \n"""
    poscar_content = []

    if '.cif' in structure_file:                            # for .cif files
        parser = CifParser(structure_file)                  # create structure from .cif using pymatgen
        structure = parser.get_structures()[0]              #
        p_content = Poscar(structure)                       # create poscar from structure
        p = str(p_content).split('\n')                      # These two lines are to remodel poscar_content to the same
        poscar_content = [string + '\n' for string in p]    # data format we get from aflow - a list of strings

    else:  # this is for aflow structure files
        with open(structure_file, 'r') as f:  # need to get POSCAR content from big structure file
            for line in f:
                if 'SPRIM' in line:  # subsection we are interested in goes after line with word "Representative"
                    count = 0
                    for line in f:  # now you are at the lines you want
                        count = count + 1
                        if 'SCONV' in line:  # Ends before line with "WYCCAR"
                            break
                        else:
                            poscar_content.append(line)
                            if count == 5:              # aflow for some reason slightly deviates from standart POSCAR
                                                        # format - the element names are missing, VASP is okay
                                                        # with this, but it is slightly inconvinient if we want to be
                                                        # consistent with other databases, so I add a dummy line
                                poscar_content.append("missing_element_names_by_Aflow\n")

    return poscar_content



