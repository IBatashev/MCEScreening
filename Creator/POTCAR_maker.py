import pandas as pd


def writer(out_path, at_type_list):
    """Creates POTCAT file by copying and joining relevant POTCARS from pots directory"""

    # pots_dir = /vol/thchem/dewijs/vasp/potcars/PBE.54/
    pots_dir = '../pots/'  # Temporary probably
    df = pd.read_csv('recommended_PAW.csv', index_col=0, sep=',')  # read the recommended potentials list

    pathlist = []
    for item in at_type_list:                      # Create a list of unique elements that are present in the compound
        POTCAR_type = df.loc[item, 'POTCAR']       # get the recommended potential type
        path = pots_dir + POTCAR_type.strip(' ') + '/POTCAR'
        pathlist.append(path)                      # create a list of paths to POTCAR files for all necessary elements

    with open(out_path+'POTCAR', 'w') as outfile:   # concatenate all the files and write out the final POTCAR
        for fname in pathlist:
            with open(fname) as infile:
                for line in infile:
                    outfile.write(line)
