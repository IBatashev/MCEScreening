import numpy as np
import pandas as pd

def writer(out_path, poscar_info):
    """Creates POTCAT file by copying and joining relevant POTCARS from pots directory"""

    # pots_dir = /vol/thchem/dewijs/vasp/potcars/PBE.54/
    pots_dir = '../pots/' # Temporary probably

    df = pd.read_csv('recommended_PAW.csv', index_col=0, sep=',')  # read the recommended potentials list

    at_num = np.fromstring(poscar_info[5], dtype=np.int, sep=' ')  # get the number of elements from POSCAR
    at_num_total = np.sum(at_num)
    at_type_list = []

    for i in range(7, 7 + int(at_num_total)):       # POSCAR files from aflowlib have 7 lines of text before atomic coordinates! so we start at line 8
        l = poscar_info[i].split()
        at_type_list.append(str(l[3]))

    pathlist = []
    for item in np.unique(at_type_list):            # Create a list of unique elements that are present in the compound
        POTCAR_type =  df.loc[item, 'POTCAR']       # get the recommended potential type
        path = pots_dir + POTCAR_type + '/POTCAR'
        pathlist.append(path)                       # create a list of paths to POTCAR files for all necessary elements

    with open(out_path+'POTCAR', 'w') as outfile:   # concaterate all the files and write out the final POTCAR
        for fname in pathlist:
            with open(fname) as infile:
                for line in infile:
                    outfile.write(line)