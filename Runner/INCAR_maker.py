import numpy as np
import pandas as pd

def writer(out_path, poscar_info):
    """Creates INCAR file based on initial POSCAR from aflow in the selected path"""

    at_num = np.fromstring(poscar_info[5], dtype=np.int, sep=' ')
    at_num_total = np.sum(at_num)
    at_type_list = []
    for i in range(7, 7 + int(at_num_total)):  # POSCAR files from aflowlib have 7 lines of text before atomic coordinates! so we start at line 8
        l = poscar_info[i].split()
        at_type_list.append(str(l[3]))


    datafile = pd.read_csv('recommended_PAW.csv', index_col=0, sep=',')
    encut_list = np.empty([0, 1])
    moment_list = np.empty([0, 1])
    for i in at_type_list:
        encut_list = np.append(encut_list, datafile.loc[i, 'ENCUT'])
        moment_list = np.append(moment_list, datafile.loc[i, 'mag_moment'])
    encut = int(1.3*np.max(encut_list)) # I set the usual "safe" increase of ENCUT by 1.3
    incar = open(out_path+'INCAR','w')
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
    "ENCUT = " + str(encut) +"\n"
    )
    incar.write(
    "MAGMOM = " + np.array_str(moment_list).strip('[').strip(']')+ "\n"
    )