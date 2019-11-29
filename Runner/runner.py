import os
import numpy as np
from ase import geometry

from Runner import INCAR_maker
from Runner import POSCAR_maker

def single_run(def_type, def_matrix):
    """Creates files for a single VASP run, executes job and cleans up"""

    path_to_calc = path + def_type + "/"
    os.mkdir(path_to_calc)
    POSCAR_maker.writer(path_to_calc, poscar_content, def_matrix)
    INCAR_maker.writer(path_to_calc, poscar_content)
    # POTCAR maker
    # job maker
    # run job
    # os.remdir after running anf copiing all info

def undeformed_lattice():
    """Creates all necessary files to run a VASP job for initial undeformed structure"""

    deformation_type = "undeformed"
    deformation_matrix = [[1,1,1],[1,1,1],[1,1,1]]
    single_run(deformation_type, deformation_matrix)

def deformed_lattice(lattice_type):
    """Takes lattice type as argument and depending on symmetry creates all necessary files to run VASP jobs for all possible independent deformations.

    For questions concerning this section one should consult bravais_lattice.pdf that contains info on lattice vectors
    for primitive cells of all 14 types. The deformation matrix is determined in a way to certainly add deformation to
    either 'a', 'b' or 'c' parameter and extra unneeded multipliers 'n' would be removed upon multiplication with actual
    lattice matrix in POSCAR_maker.py. e.g.:
     deformation matrix          primitive lattice matrix                                       deformed lattice matrix
    [[n,n,n],[1,1,1],[1,1,1]] * [[3,0,3],[2,0,2],[0,0,1]] = [[3*n,0*n,3*n],[2,0,2],[0,0,1]] = [[3n,0,3n],[2,0,2],[0,0,1]]
    """

    # "CUB"     1 Simple cubic
    # "FCC"     2 Face centered cubic
    # "BCC"     3 Body centered cubic
    # "HEX"     4 Hexagonal close packed
    # "TET"     5 Simple tetragonal
    # "BCT"     6 Body centered tetragonal
    # "RHL"     7 Rhombohedral
    # "ORC"     8 Simple orthogonal
    # "ORCC"    9 Base centered orthorhombic
    # "ORCI"    10 Body centered orthorhombic
    # "ORCF"    11 Face centered orthorhombic
    # "MCL"     12 Simple monoclinic
    # "MCLC"    13 Base centered monoclinic
    # "TRI"     14 Triclinic

    if lattice_type == "CUB" or lattice_type == "FCC" or lattice_type == "BCC":
        deformation_types = 'a'
        for deformation_type in deformation_types:
            if deformation_type == "a":
                deformation_matrix = np.array([[n,n,n],[n, n,n],[n,n,n]])
                single_run(deformation_type, deformation_matrix)

    elif lattice_type == "HEX":
        deformation_types = 'ac'
        for deformation_type in deformation_types:
            if deformation_type == "a":
                deformation_matrix = np.array([[n,n,0],[n,n,0],[0,0,1]])
                single_run(deformation_type, deformation_matrix)
            if deformation_type == "c":
                deformation_matrix = np.array([[1,1,0],[1,1,0],[0,0,n]])
                single_run(deformation_type, deformation_matrix)

    elif lattice_type == "RHL":                                                             # very confused by rhombohedral lattice
        print("very confused by RHL need to think and write this part")

    elif lattice_type == "TET" or lattice_type == "BCT":
        deformation_types = 'ac'
        for deformation_type in deformation_types:
            if deformation_type == "a":
                deformation_matrix = np.array([[n,n,1],[n,n,1],[n,n,1]])
                single_run(deformation_type, deformation_matrix)
            if deformation_type == "c":
                deformation_matrix = np.array([[1,1,n],[1,1,n],[1,1,n]])
                single_run(deformation_type, deformation_matrix)

    elif lattice_type == "ORC" or lattice_type == "ORCC" or lattice_type == "ORCI" or lattice_type == "ORCF":
        # the deformation matrixies I now use should properly cover all 4 types of primitive orthorhombic lattices
        deformation_types = 'abc'
        for deformation_type in deformation_types:
            if deformation_type == "a":
                deformation_matrix = np.array([[n,1,1],[n,1,1],[n,1,1]])
                single_run(deformation_type, deformation_matrix)
            if deformation_type == "b":
                deformation_matrix = np.array([[1,n,1],[1,n,1],[1,n,1]])
                single_run(deformation_type, deformation_matrix)
            if deformation_type == "c":
                deformation_matrix = np.array([[1,1,n],[1,1,n],[1,1,n]])
                single_run(deformation_type, deformation_matrix)

    elif lattice_type == "MCL" or lattice_type == "MCLC":
        deformation_types = 'abc'
        for deformation_type in deformation_types:
            if deformation_type == "a":
                deformation_matrix = np.array([[n,1,1],[n,1,1],[1,1,1]])
                single_run(deformation_type, deformation_matrix)
            if deformation_type == "b":
                deformation_matrix = np.array([[1,n,1],[1,n,1],[1,1,1]])
                single_run(deformation_type, deformation_matrix)
            if deformation_type == "c":
                deformation_matrix = np.array([[1,1,1],[1,1,1],[n,n,n]]) # Aflow for some reason is different for a3 = ccosβx^+csinβz so 'a1' and 'a2' are switched leading to 'a' and 'b' being switched
                single_run(deformation_type, deformation_matrix)

    elif lattice_type == "TRI":
        deformation_types = 'abc'
        for deformation_type in deformation_types:
            if deformation_type == "a":
                deformation_matrix = np.array([[n,0,0],[1,1,0],[1,1,1]])
                single_run(deformation_type, deformation_matrix)
            if deformation_type == "b":
                deformation_matrix = np.array([[1,0,0],[n,n,0],[1,1,1]])
                single_run(deformation_type, deformation_matrix)
            if deformation_type == "c":
                deformation_matrix = np.array([[1,0,0],[1,1,0],[n,n,n]])
                single_run(deformation_type, deformation_matrix)

    else: print("Error! Unknown Lattice type for "+str(ID))

ID = 1770           # at the begining script gets ID of the structure we are going to work with as input argument
path = '../test/'   # folder where all subfolders for a single ID will be created/executed/cleaned - working directory
n = 1.2             # deformation coefficient
poscar_file = '../Database/datadir/'+str(ID)

### Take all information from poscar as list of strings:
poscar = open(poscar_file, "r")
poscar_content = []
for line in poscar:
    poscar_content.append(line)
poscar.close()

### Use ase to get Bravais lattice type from lattice in poscar
lattice_matrix = np.zeros([3, 3])
for i in range(2, 5):
    lattice_matrix[i - 2, :] = np.fromstring(poscar_content[i], dtype=np.float, sep=' ')
read_lattice_matrix = geometry.Cell.new(lattice_matrix)
lat_type = str(read_lattice_matrix.get_bravais_lattice()).split("(")[0]

### Create a folder for a job without any deformation and run it
undeformed_lattice()

### Create the proper number of folders for all possible deformations according to Bravais lattice type and run them
deformed_lattice(lat_type)