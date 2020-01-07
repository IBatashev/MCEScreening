import pandas as pd
import os
import ast

# Made this since I accidentially downloaded some compounds with U, Th and Pu and had no need to keep the entries for them...

def delete(datalist, datadir):
    """Completely removes entry from datalist and deletes all corresponding files from datadir
    WARNING use with care - unreversable, so better have backup of original datafiles somewhere..."""

    element = 'Pu'

    df = pd.read_csv(datalist, index_col=0, sep=',')
    for item in df.index.tolist():
        elem_list = str(df.loc[item, 'species'])
        if element in elem_list:
            os.remove(datadir +'_incars/' + str(item))
            os.remove(datadir +'_aflow/' + str(item))
            os.remove(datadir +'_structure/' + str(item))
            os.remove(datadir +'_structure_relaxed/' + str(item))
            df = df.drop([item], axis=0)
    df.to_csv(datalist.replace(".csv",'_cleaned.csv'))

delete('../Database/datalist.csv','../Database/datadir')