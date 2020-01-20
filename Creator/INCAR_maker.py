import numpy as np
import pandas as pd


def writer(out_path, poscar_info, at_type_list):
    """Creates INCAR file based on initial POSCAR from aflow in the selected path"""

    at_num = np.fromstring(poscar_info[5], dtype=np.int, sep=' ')
    datafile = pd.read_csv('recommended_PAW.csv', index_col=0, sep=',') # csv file containing all ENCUTS (ENMAX from POTCARs) and moments (from "intuition")
    encut_list = np.empty([0, 1])
    magmom = ''

    for count, elem in enumerate(at_type_list):
        encut_list = np.append(encut_list, datafile.loc[elem, 'ENCUT'])
        moment = datafile.loc[elem, 'MOMENT']
        magmom = magmom + (str(at_num[count]) + '*' + str(moment)+' ')  # make a string of number of atoms multiplied by their moments
    encut = int(1.3 * np.max(encut_list))                               # biggest ENCUT among all elements in the composition with the usual "safe" increase of 1.3

    incar = open(out_path + 'INCAR', 'w')
    incar.write(
        "ALGO = Fast \n"
        "PREC = Accurate \n"
        "ISIF = 1 \n"
        "ISPIN = 2 \n"
        "ISMEAR = -5 \n"
        "NELM = 100 \n"
        "SIGMA = 0.05 \n"
        "LORBIT = 10 \n"
        "IBRION = 2 \n"
        "NSW = 100 \n"
        "EDIFF = 0.00001 \n"
        "NCORE = 4 \n"
        "LCHARG = .FALSE. \n"
        "LWAVE = .FALSE. \n"
    )
    incar.write(
        "ENCUT = " + str(encut) + "\n"
    )
    incar.write(
        "MAGMOM = " + magmom
    )