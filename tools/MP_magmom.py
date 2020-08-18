import numpy as np
import os
import pandas as pd
from pymatgen.io.vasp import Poscar
from Tools import POSCAR_reader


def moment_parser(item):
    """Reads moments from MAGMOM tag from datalist"""

    magmom_string = (df.loc[item, 'magmom']).strip('[').strip(']').strip(' ')
    magmoms_list = magmom_string.split(';')
    mom = np.array(magmoms_list, float)
    return mom


def element_parser(item):
    """Returns list of all elemets from structure file"""

    structure_file = wdatadir + str(item) + '.cif'
    poscar_content = POSCAR_reader.read(structure_file)
    poscar_string = ''.join(poscar_content)  # merging poscar content in a single string so pymatgen can read it
    poscar = Poscar.from_string(poscar_string)  # using pymatgen to acquire our structure from poscar content
    structure = poscar.structure
    elem_string = str(structure.species).strip('[').strip(']').replace("Element", "").replace(" ", "")
    elem_list = elem_string.split(',')
    elem = np.array(elem_list)
    return elem


def write(ID):
    """Appends magnetic moment to the corresponding element column in csv file """

    momlist = moment_parser(ID)
    elemlist = element_parser(ID)

    # el_df = pd.read_csv('elem_moments.csv', sep=',',index_col=0)

    for num, i in enumerate(elemlist, 0):

    ### works but only provide info for average
    #     el_df.loc[i, 'moment_sum'] = (momlist[num])+el_df.loc[i,'moment_sum']
    #     el_df.loc[i, 'moment_count'] = 1+el_df.loc[i,'moment_count']
    # el_df.to_csv("elem_moments.csv")
    ### create shitload of txts
        with open('./moments_by_elements/' + i + '.txt', 'a') as f:
            f.write(str(momlist[num]) + '\n')


def averager(dir):
    with open('avg_moments.csv', 'a') as f:  # we create a csv file to write info into
        f.write("element,average_moment\n")

    filelist = os.listdir(dir)
    for i in filelist:
        element = i.strip('.txt')
        moments = np.loadtxt(dir + i)
        avg_mom = np.average(moments)

        newrow = str(element + ',' + str(avg_mom) + '\n')
        with open('avg_moments.csv', 'a', encoding="utf-8") as f:
            f.write(newrow)



wdatadir = '../Database/MP/datadir/'
wdatalist = '../Database/MP/datalist_lattfix_updated_sieved.mag.field_sieved.mag.sites_no.duplicates.csv'
directory = './moments_by_elements/'

df = pd.read_csv(wdatalist, index_col=0, sep=',')

for item in df.index.tolist():
    mom_arr = moment_parser(item)
    total = mom_arr.sum()
    df.loc[item, 'sum_of_moments'] = total
df.to_csv(wdatalist.replace(".csv", '_mom_summed' + '.csv'))

#     try:
#         write(item)
#     except:
#         with open('./faillist.txt', 'a') as f:
#             f.write(str(item)+'\n')
#

# averager(directory)

