# v1.0 11-10-2019
# python=3
import numpy as np
import ase
import ase.cell as cell

def KSTR_maker(ID):
    """This makes a KSTR file from POSCAR.
    Technically works, but not tested, as EMTO calculations are not a priority for now.
    So this script could definetly be optimised."""

    lat_list = [
        "CUB",                        # 1 Simple cubic
        "FCC",                        # 2 Face centered cubic
        "BCC",                        # 3 Body centered cubic
        "HEX",                        # 4 Hexagonal close packed
        "TET",                        # 5 Simple tetragonal
        "BCT",                        # 6 Body centered tetragonal
        "RHL",                        # 7 Rhombohedral
        "ORC",                        # 8 Simple orthogonal
        "ORCC",                       # 9 Base centered orthorhombic
        "ORCI",                       # 10 Body centered orthorhombic
        "ORCF",                       # 11 Face centered orthorhombic
        "MCL",                        # 12 Simple monoclinic
        "MCLC",                       # 13 Base centered monoclinic
        "TRI"                         # 14 Triclinic
    ]

    ### Extracting data from POSCAR ###
    poscar = open("./data/"+ID, "r")    # opening file from path
    lines = []
    for line in poscar:                 # reading all text from POSCAR file
        lines.append(line)
    poscar.close()

    # Getting lattice matrix from POSCAR
    lat_mat = np.zeros([3, 3])
    for i in range(2, 5):
        lat_mat[i-2, :] = np.fromstring(lines[i], dtype=np.float, sep=' ')

    # Getting Bravias lattice type using ase and matrix from poscar
    lat_cell= ase.cell.Cell.new(lat_mat)
    lat_type_name = str(ase.cell.Cell.get_bravais_lattice(lat_cell)).split('(')[0]
    lat_type = lat_list.index(lat_type_name)+1

    # Getting total number of atoms from POSCAR
    num_at = sum(list((map(int, lines[5].split()))))

    # Getting all coordinates
    coord = np.zeros([num_at, 3])
    i =0
    while i < num_at:
        l = lines[7+i].split()
        coord[i, 0] = l[0]
        coord[i, 1] = l[1]
        coord[i, 2] = l[2]
        i = i+1

    # WRITING KSTR file
    output_path = "./"+ID+".kstr"
    KSTR = open(output_path, 'w')

    KSTR.write(
        "KSTR HP......=N             09 Mar 99\n"
        "JOBNAM...="+ID+" MSGL.= 1 MODE...=B STORE..=Y HIGH...=Y\n"
        "FOR001=smx/\n"
        "FOR006=./\n"
        "Slope matrices , can_i_write_anything_here?\n"
        "NL.....= 4 NLH...= 9 NLW...= 7 NDER..= 6 ITRANS= 3 NPRN..= 0\n"
        "(K*W)^2..= 0.000000 DMAX....= 2.5000 RWATS...= 0.10\n" # Carefull with DMAX, in the end we still need to make it automatic
        "NQ3...= "+str(num_at)+" LAT...= "+str(lat_type)+" IPRIM.= 0 NGHBP.=13 NQR2..= 0\n" # NQ3  - this is a number of atoms in total, right?
        "A........=1.0 B.......=1.0 C.......=1.0\n"
    )
    KSTR.write(
        "BSX......=%.8f BSY.....=%.8f BSZ.....=%.8f\n"  # Does precision matter for EMTO? POSCAR from AFLOW gives 14 digits after decimal point
        "BSX......=%.8f BSY.....=%.8f BSZ.....=%.8f\n"  # I set here up to 8th, but obviously can include more digits
        "BSX......=%.8f BSY.....=%.8f BSZ.....=%.8f\n"
        %(lat_mat[0, 0], lat_mat[0, 1], lat_mat[0, 2],
           lat_mat[1, 0], lat_mat[1, 1], lat_mat[1, 2],
           lat_mat[2, 0], lat_mat[2, 1], lat_mat[2, 2])
    )
    i=0
    while i < num_at:
        KSTR.write(
            "QX.......=%.8f QY......=%.8f QZ......=%.8f\n" %(coord[i,0],coord[i,1],coord[i,2])
        )
        i = i+ 1
    KSTR.write(
        "a/w(.1)..= 0.70 0.70 0.70 0.70\n"
        "a/w(.2)..= 0.70 0.70 0.70 0.70\n"
        "a/w(.3)..= 0.70 0.70 0.70 0.70\n"
        "a/w(.4)..= 0.70 0.70 0.70 0.70\n"
        "NL_mdl.= 9\n"
        "LAMDA....= 2.5000 AMAX....= 4.5000 BMAX....= 4.5000"
    )