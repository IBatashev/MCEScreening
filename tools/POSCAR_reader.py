from pymatgen.io.cif import CifParser
from pymatgen.io.vasp import Poscar
from pymatgen.core import Structure
import warnings
import re

warnings.filterwarnings("ignore")  # slightly dirty, but simple way to disable pymatgen complaining about lack
                                   # of element labels in POSCARs from aflowlib "UserWarning: Elements in POSCAR cannot be determined."
                                   # Warning ! This disables _ALL_ warnings.


def read(structure_file):
    """Takes structure file as input and returns section containing POSCAR as list of strings ending with \n"""
    poscar_content = []

    if '.cif' in structure_file:                            # for .cif files
        #parser = CifParser(structure_file)                  # create structure from .cif using pymatgen
        #structure = parser.get_structures()[0]              #
        mystructure = Structure.from_file(structure_file)
        # trying to deal with partial occupancies:

        if mystructure.is_ordered:
            p_content = Poscar(mystructure)                       # create poscar from structure
            p = str(p_content).split('\n')                      # These two lines are to remodel poscar_content to the same
            poscar_content = [string + '\n' for string in p]    # data format we get from aflow - a list of strings
        else:
            print('disordered structure')


    else:
        with open(structure_file) as f:
            first_line = f.readline()
            if 'VERSION' in first_line:

                # this is for aflow structure files
                with open(structure_file, 'r') as f:  # need to get POSCAR content from big structure file
                    for line in f:
                        if 'WYCCAR' in line:
                            compound = (next(f).split()[0].replace(".", ""))
                            element_string = re.sub(r'\d+', '', compound)
                            output = str((re.findall('[A-Z][^A-Z]*', element_string))).strip('[').strip(']').replace("'", "").replace(",", "")
                        if 'SPRIM' in line:
                            count = 0
                            for line in f:
                                count = count + 1
                                if 'SCONV' in line:
                                    break
                                else:
                                    poscar_content.append(line)
                                    if count == 5:              # aflow for some reason slightly deviates from standart POSCAR
                                                                # format - the element names are missing, so we get them from other part of structure file
                                        poscar_content.append(output + "\n")
            else:
                with open(structure_file, 'r') as f:  # need to get POSCAR content from big structure file
                    for line in f:
                        poscar_content.append(line)

    return poscar_content

