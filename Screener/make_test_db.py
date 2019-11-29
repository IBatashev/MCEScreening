import pandas as pd

def Make_Sample_set_for_MCEScreening_tests():
    # make this shit automatic also make it move files to testDatabase folder
    df = pd.read_csv('Alldata.csv', index_col=0, sep=',')
    monoclinic = df[df['lattice_system'].str.match('monoclinic')]       # no preset
    triclinic = df[df['lattice_system'].str.match('triclinic')]         # no preset
    rhombohedral = df[df['lattice_system'].str.match('rhombohedral')]   # no preset
    hexagonal = df[df['lattice_system'].str.match('hexagonal')]         # Fe2P
    orthorhombic = df[df['lattice_system'].str.match('orthorhombic')]   # MnCoP
    cubic = df[df['lattice_system'].str.match('cubic')]                 # FeRh
    tetragonal = df[df['lattice_system'].str.match('tetragonal')]       # LaFe9Si4

    k = tetragonal
    out = k.sample(n=4, random_state=1)
    out.to_csv('100/tetragonal.csv')
