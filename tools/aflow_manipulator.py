import numpy as np
import linecache
import pandas as pd

# Various functions to work with aflow datalist and datadirs.


def INCAR_parser(incar_file):
    """Reads moments from MAGMOM tag from INCAR files"""

    incar = open(incar_file, "r")
    for line in incar:
        if 'MAGMOM=' in line:
            m = str.split(line, '= ')
            m = m[1].split('#')
            magmom = m[0].split()
    #     # else: magmom = ''
    incar.close()
    num = []
    mom = []

    for elem in magmom:
        if '*' in elem:
            num = np.append(num, elem.split('*')[0])
            mom = np.append(mom, elem.split("*")[1])
        else:
            mom = np.append(mom, elem)
            num = np.append(num, '1')
        num = np.ndarray.astype(num, int)
        mom = np.ndarray.astype(mom, float)
    magnetic = []
    for i in range(0, len(num)):
        for k in range(0, num[i]):
            magnetic = np.append(magnetic, mom[i])
    return magnetic
    # if no magmom return exception


def element_parser(structure_file):
    """Returns list of all elemets from structure file"""
    structure = open(structure_file, "r")
    ii = 0
    for num, line in enumerate(structure, 1):
        if 'Direct' in line:
            if ii == 1:
                linecount = int((line.split('(')[1].split(')'))[0])
                break
            else:ii = ii+1
    i = 0
    elemlist = []
    while i < linecount:
        i = i+1
        elemline = linecache.getline(structure_file, num+i)
        elem = elemline.split()[3]
        elemlist = np.append(elemlist, elem)
    return elemlist


def write(ID):
    """Appends magnetic moment to the corresponding element column in csv file """
    momlist = (INCAR_parser("../Database/datadir_incars/"+str(ID)))
    elemlist = (element_parser("../Database/datadir_structure/"+str(ID)))
    el_df = pd.read_csv('elem_moments.csv', sep=',',index_col=0)
    for num, i in enumerate(elemlist,0):
    ### works but only provide info for average
    #     el_df.loc[i, 'moment_sum'] = (momlist[num])+el_df.loc[i,'moment_sum']
    #     el_df.loc[i, 'moment_count'] = 1+el_df.loc[i,'moment_count']
    # el_df.to_csv("elem_moments.csv")
    ### create shitload of txts
        with open('./moments_by_elements/' + i + '.txt', 'a') as f:
            f.write(str(momlist[num]) + '\n')


df = pd.read_csv("../Database/datalist.csv", index_col=0, sep=',')
for item in df.index.tolist():
    try:
        write(item)
        print(item)
    except:
        with open('./faillist.txt', 'a') as f:
            f.write(str(item)+'\n')


# print(INCAR_parser('../Database/datadir_incars/1'))
