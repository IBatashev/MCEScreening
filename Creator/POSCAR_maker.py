from pymatgen.io.vasp import Poscar
import numpy as np


def writer(out_path, poscar_content, deformation):
    """Creates POSCAR in the selected path using initial
    POSCAR content from aflow structure file and chosen deformation"""

    poscar_string = ''.join(poscar_content)       # merging poscar content in a single string so pymatgen can read it
    poscar = Poscar.from_string(poscar_string)    # using pymatgen to acquire our structure from poscar content
    structure = poscar.structure

    lattice_string = str(structure.lattice).replace('\n', ' ')          # convert pymatgen lattice to string
    lattice = np.fromstring(lattice_string, sep=' ').reshape(3, 3)      # convert resulting string to numpy array
    lattice_deformed = np.multiply(lattice, deformation)                # apply deformation
    structure.lattice.__init__(lattice_deformed)                        # update pymatgen structure with new lattice

    filename = 'POSCAR'
    o = open(out_path+filename, "w+")
    o.write(poscar_content[0])                    # Writing poscar header
    o.write('1.0\n')                              # Writing constant by which all lattice parameters are multiplied - hopefully pymatgen will always give lattice matrix without any scaling so we can just put 1.0 here
    o.write(str(structure.lattice) + '\n')        # writing lattice matrix from pymatgen
    o.write(poscar_content[5])                    # writing type of atoms
    o.write(poscar_content[6])                    # writing number of atoms - copied from original poscar
    o.write('Direct' + '\n')                      # writing type of coordinates, again we will always get direct coords from pymatgen, so we must have 'Direct' here
    o.write(' ' + str(structure.frac_coords).replace('[', '').replace(']', ''))  # writing direct coordinates from pymatgen


def writer_MP(out_path, poscar_content, lattice):
    """Creates POSCAR in the selected path using initial
    POSCAR content from MP structure file and chosen deformation"""

    poscar_string = ''.join(poscar_content)       # merging poscar content in a single string so pymatgen can read it
    poscar = Poscar.from_string(poscar_string)    # using pymatgen to acquire our structure from poscar content
    structure = poscar.structure

    lattice_string = str(lattice).replace('\n', ' ')  # convert pymatgen lattice to string
    update_lattice = np.fromstring(lattice_string, sep=' ').reshape(3, 3)  # convert resulting string to numpy array

    structure.lattice.__init__(update_lattice)                        # update pymatgen structure with new lattice

    filename = 'POSCAR'
    o = open(out_path+filename, "w+")
    o.write(poscar_content[0])                    # Writing poscar header
    o.write('1.0\n')                              # Writing constant by which all lattice parameters are multiplied - hopefully pymatgen will always give lattice matrix without any scaling so we can just put 1.0 here
    o.write(str(structure.lattice) + '\n')        # writing lattice matrix from pymatgen
    o.write(poscar_content[5])                    # writing type of atoms
    o.write(poscar_content[6])                    # writing number of atoms - copied from original poscar
    o.write('Direct' + '\n')                      # writing type of coordinates, again we will always get direct coords from pymatgen, so we must have 'Direct' here
    o.write(' ' + str(structure.frac_coords).replace('[', '').replace(']', ''))  # writing direct coordinates from pymatgen

