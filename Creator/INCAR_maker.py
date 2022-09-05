import numpy as np
import pandas as pd
import math


def writer(out_path, poscar_info, at_type_list, calculation_type, moments, EDIFF_scaling=False):
    """Creates INCAR file based on initial POSCAR from aflow in the selected path"""

    at_num = np.fromstring(poscar_info[6], dtype=np.int, sep=' ')
    datafile = pd.read_csv('recommended_PAW.csv', index_col=0, sep=',') # csv file containing all ENCUTS (ENMAX from POTCARs) and moments (from "intuition")
    encut_list = np.empty([0, 1])
    magmom = ''
    number_of_atoms = at_num.sum()
    EDIFF = (math.sqrt(number_of_atoms)*0.0001)

    if moments == '':
        for count, elem in enumerate(at_type_list):
            encut_list = np.append(encut_list, datafile.loc[elem, 'ENCUT'])
            moment = datafile.loc[elem, 'MOMENT']
            magmom = magmom + (str(at_num[count]) + '*' + str(moment)+' ')  # make a string of number of atoms multiplied by their moments
    # else:
    #     for count, elem in enumerate(at_type_list):
    #         encut_list = np.append(encut_list, datafile.loc[elem, 'ENCUT'])
    #     magmom = moments

    encut = int(1.3 * np.max(encut_list))                               # biggest ENCUT among all elements in the composition with the usual "safe" increase of 1.3

    incar = open(out_path + 'INCAR', 'w')
    incar.write(
        "ALGO = Fast \n"
        "PREC = Accurate \n"
        "ISIF = 2 \n"
        "ISPIN = 2 \n"
        "ISMEAR = 2  \n"
        "NELM = 50 \n"
        "SIGMA = 0.2 \n"
        "KSPACING = 0.5 \n"
        "LORBIT = 10 \n"
        "NCORE = 4 \n"
        "LASPH = .TRUE. \n"
        "GGA_COMPAT = .FALSE. \n"
        "LREAL = A \n"

        # "AMIX = 0.1 \n"
        # "BMIX = 0.00001 \n"
        # "AMIX_MAG = 0.2 \n"
        # "BMIX_MAG = 0.00001 \n"
    )
    incar.write("ENCUT = " + str(encut) + " \n")
    incar.write("MAGMOM = " + magmom + " \n")
    if calculation_type == 'Applied_Field':
        incar.write("BEXT = -0.05 \n")
    if EDIFF_scaling == True:
        incar.write("EDIFF = " + str(EDIFF) + " \n")
    else:
        incar.write("EDIFF = 0.00001 \n")


def writer_second(out_path, poscar_info, at_type_list, deformation):
    """Creates INCAR file based on initial POSCAR from aflow in the selected path"""

    at_num = np.fromstring(poscar_info[6], dtype=np.int, sep=' ')
    datafile = pd.read_csv('recommended_PAW.csv', index_col=0, sep=',') # csv file containing all ENCUTS (ENMAX from POTCARs) and moments (from "intuition")
    encut_list = np.empty([0, 1])
    magmom = ''

    for count, elem in enumerate(at_type_list):
        encut_list = np.append(encut_list, datafile.loc[elem, 'ENCUT'])
        moment = datafile.loc[elem, 'MOMENT']
        magmom = magmom + (str(at_num[count]) + '*' + str(moment)+' ')  # make a string of number of atoms multiplied by their moments
    encut = int(1.5 * np.max(encut_list))                               # biggest ENCUT among all elements in the composition with a big increase of 1.5

    incar = open(out_path + 'INCAR', 'w')
    incar.write(
        "ALGO = Fast \n"
        "PREC = Accurate \n"
        "ISIF = 3 \n"
        "ISPIN = 2 \n"
        "ISMEAR = 2  \n"
        "NELM = 100 \n"
        "SIGMA = 0.2 \n"
        "KSPACING = 0.5 \n"
        "LORBIT = 10 \n"
        "EDIFF = 0.00001 \n"
        "NCORE = 4 \n"
        "LASPH = .TRUE. \n"
        "GGA_COMPAT = .FALSE. \n"
        #"LCHARG = .FALSE. \n"
        #"LWAVE = .FALSE. \n"
        "LREAL = A \n"
        "NSW = 100 \n"
        "IBRION = 2 \n"
    )
    incar.write(
        "ENCUT = " + str(encut) + " \n"
    )
    incar.write(
        "MAGMOM = " + magmom + " \n"
    )
    if deformation == 'BEXT':
        incar.write(
            "BEXT = -0.05 \n"  # uncomment for applied filed
        )
