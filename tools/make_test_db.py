import pandas as pd
import shutil
import os

def Make_Sample_set_for_MCEScreening_tests(random_seed=1, sample_size=4):
    """"Takes a set number (sample_size) of random compounds for each lattice type from full database and creates a new
    datalist and datadir to use for quick testing of screening approach"""

    new_datadir =  '../Database/TestDB/'
    lattice_systems = ['monoclinic','triclinic','rhombohedral','hexagonal','orthorhombic','cubic','tetragonal'] # Here we choose which lattice sistems would be selected from full database
    df = pd.read_csv('../Database/datalist_updated_sieved.mag.field_sieved.mag.sites_no.duplicates.csv', index_col=0, sep=',')
    df_out = pd.DataFrame()                                                                     # Prepare new empty pandas dataframe
    if os.path.exists(new_datadir):                                                             # Prepare new datadir folder
        shutil.rmtree(new_datadir)                                                              # (cleans old one if it already existed)
    os.makedirs(new_datadir)
    for lattice in lattice_systems:                                                             # For each lattice system choose a number of random compounds
        lat_type_subset = df[df['lattice_system'].str.match(lattice)]
        random_sample = lat_type_subset.sample(n=sample_size, random_state=random_seed)
        df_out = pd.concat([df_out, random_sample])                                             # Fill the new dataframe with resulting sample
    df_out.to_csv('../Database/datalist_TestDB.csv')                                            # Write a new .csv file
    for item in df_out.index.tolist():                                                          # Iterate over IDs in our new small datalist
        shutil.copy('../Database/datadir_structure_relaxed/'+str(item),new_datadir+str(item))                     # and copy POSCAR files from full datadir to new smaller one

Make_Sample_set_for_MCEScreening_tests(1, 4)

# Fe2P 6295
# FeRh
# LaFeSi 5565
# LaSrMnO3
# MnCoP