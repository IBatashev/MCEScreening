from pymatgen.io.vasp import Poscar


def writer(out_path, poscar_content, deformation):
    """Creates POSCAR in the selected path using initial POSCAR content from aflow structure file and pymatgen deformation"""

    poscar_string = ''.join(poscar_content)       # merging poscar content in a single string so pymatgen can read it
    poscar = Poscar.from_string(poscar_string)    # using pymatgen to acquire our structure from poscar content
    structure = poscar.structure
    structure.apply_strain(deformation)           # adding deformation to the lattice with pymatgen

    filename = 'POSCAR'
    o = open(out_path+filename, "w+")
    o.write(poscar_content[0])                    # Writing poscar header
    o.write('1.0')                                # Writing constant by which all lattice parameters are multiplied - hopefully pymatgen will always give lattice matrix without any scaling so we can just put 1.0 here
    o.write(str(structure.lattice) + '\n')        # writing lattice matrix from pymatgen
    o.write(poscar_content[5])                    # writing number of atoms - copied from original poscar
    o.write('Direct' + '\n')                             # writing type of coordinates, again we will always get direct coords from pymatgen, so we must have 'Direct' here
    o.write(' ' + str(structure.frac_coords).replace('[', '').replace(']', ''))  # writing get direct coordinates from pymatgen
