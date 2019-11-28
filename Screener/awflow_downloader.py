from aflow import *
import time

result = search(catalog='icsd', batch_size=100
    ).filter( (K.species == "Mn") |
              (K.species == "Fe") |
              (K.species == "Co") |
              (K.species == "Ni") |
              (K.species == "Cu")
    ).select(
             ~
              (K.species == "Rb") &~
              (K.species == "Cs") &~
              (K.species == "Ba") &~
              (K.species == "Re") &~
              (K.species == "Os") &~
              (K.species == "Ir") &~
              (K.species == "Pt") &~
              (K.species == "Au") &~
              (K.species == "Hg") &~
              (K.species == "Tl") &~
              (K.species == "Pb") &~
              (K.species == "Bi") &~
              (K.species == "Tb") &~
              (K.species == "Dy") &~
              (K.species == "Ho") &~
              (K.species == "Er") &~
              (K.species == "Tm") &~
              (K.species == "Yb") &~
              (K.species == "Lu") &~
              (K.species == "As") &~
              (K.species == "He") &~
              (K.species == "Ne") &~
              (K.species == "Ar") &~
              (K.species == "Kr") &~
              (K.species == "Xe")
    ).select(K.files == "CONTCAR.relax.vasp"
             # This is only to look at entries that have POSCAR.
             # This selection is not really necessary as it seems that all entries in FULL aflolib have POSCAR.relax structure file
             # (165705 hits with; same 165705 hits without)
    ).select((K.spin_cell > 0) | (K.spin_cell < 0)
             # only find entries with non-zero moment not sure if one need to consider <0 though
    ).select(K.enthalpy_cell < 0)
             # entalpy must me negative - othervise structure would not form in real life


totalN = len(result)
print(totalN)

def downloader(counter=0, default_decounter=250):
    """Function used to download poscar files and fill up datalist with information from Aflow
    Works with search result, has two parameters:
    counter - The search entry where we start downloading from (technically we always want zero, only added for flexebility)
    decounter - number of entries we download in a batch before waiting for some time in order not to overburden aflow with requests"""

    with open('./datalist.csv', 'a') as f:        # we create a csv file to write info into later
        f.write("ID,compound,lattice_system_relax,spacegroup_relax,volume_cell,spin_cell,mag_field,mag_sites")
    decounter = default_decounter
    while counter <= totalN:
        result[counter].files["CONTCAR.relax.vasp"]("./data/"+str(counter))
        newrow = str(counter)+','+ str(result[counter].compound)+','+str(result[counter].lattice_system_relax.strip('\n'))+','+ str(result[counter].spacegroup_relax)+','+str(result[counter].volume_cell)+','+ str(result[counter].spin_cell)+','+'0,0\n'
        with open('./datalist.csv', 'a') as f:
            f.write(newrow)
        if decounter == 0:
            print(counter, '   ', (counter/totalN)*100, '% done') # just a simple progress indicator (may try some fancy stuff later...)
            time.sleep(300)
            decounter = default_decounter
        counter = counter + 1
        decounter = decounter - 1
