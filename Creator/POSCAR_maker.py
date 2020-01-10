import numpy as np

def writer(out_path, poscar_info, deformation_matrix):
    """Creates POSCAR using initial POSCAR from aflow and deformation matrix in the selected path"""

    lattice_matrix = np.zeros([3, 3])
    for i in range(2, 5):
        lattice_matrix[i - 2, :] = np.fromstring(poscar_info[i], dtype=np.float, sep=' ')
    lattice_matrix = np.multiply(lattice_matrix, deformation_matrix) # adding deformation to the lattice
    filename = 'POSCAR'
    o = open(out_path+filename, "w+")
    o.write(poscar_info[0])                     # Writing poscar header
    o.write(poscar_info[1])                     # Writing constant by which all lattice parameters are multiplied sometimes
                                                # in aflow it doesn't equal 1.0 - yet it doesn't seem to indicate a deformation,
                                                # just a stupid way the database keeps files in.
                                                # For now I mark them as warning in screener.py - but maybe I should be remaking POSCARs
                                                # so that they are consistent with this constang = 1.0 and changing lattice matrix below accordingly
    o.write('   %.14f  %.14f  %.14f\n' % ((lattice_matrix[0,0]), (lattice_matrix[0,1]), (lattice_matrix[0,2])))
    o.write('   %.14f  %.14f  %.14f\n' % ((lattice_matrix[1,0]), (lattice_matrix[1,1]), (lattice_matrix[1,2])))
    o.write('   %.14f  %.14f  %.14f\n' % ((lattice_matrix[2,0]), (lattice_matrix[2,1]), (lattice_matrix[2,2])))
    o.write(poscar_info[5])                     # writing number of atoms - copied from original poscar
    o.write(poscar_info[6])                     # writing type of coordinates -Perhaps a check here is necessary to make sure they are direct - copied from original poscar
    at_num = np.fromstring(poscar_info[5], dtype=np.int, sep=' ')
    at_num_total = np.sum(at_num)
    for i in range(7, 7+int(at_num_total), 1):  # writing all atomic coordinates - copied from original poscar
        o.write(poscar_info[i])