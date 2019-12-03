from aflow import *
import time

result = search(catalog='icsd', batch_size=100
    ).filter(
             ### Must contain at least on of these:
             ((K.species == "Mn") |
              (K.species == "Fe") |
              (K.species == "Co") |
              (K.species == "Ni") |
              (K.species == "Cu"))
             &
             ### Does not contain any of these because:
             ### Expensive or limited in supply:
             (K.species != "Re") &
             (K.species != "Os") &
             (K.species != "Ir") &
             (K.species != "Pt") &
             (K.species != "Au") &
             (K.species != "In") &
             ### Health Hazard:
             (K.species != "Be") &
             (K.species != "As") &
             (K.species != "Cd") &
             (K.species != "Cs") &
             (K.species != "Ba") &
             (K.species != "Hg") &
             (K.species != "Tl") &
             (K.species != "Pb") &
             ### Noble gases:
             # (K.species != "He") & only reason He is not sorted out is because i reached maximum possible querry criteria,
             # and there aren't any compounds with He anyway so it is safe to leave and include one more important instead
             (K.species != "Ne") &
             (K.species != "Ar") &
             (K.species != "Kr") &
             (K.species != "Xe")
             #  Elements afte Bi not mentioned this list are not in aflow database at all.
             # It contains only first 83 elements, so they are sorted out by default.
      ).filter(K.files == "CONTCAR.relax.vasp"
             # This is only to look at entries that have POSCAR.
             # This selection is not really necessary as it seems that all entries in FULL aflolib have POSCAR.relax structure file
             # (165705 hits with; same 165705 hits without)
      ).filter((K.spin_cell > 0) | (K.spin_cell < 0)
    #          # only find entries with non-zero moment not sure if one need to consider <0 though
      ).filter(K.enthalpy_cell < 0)
             # entalpy must me negative - othervise structure would not form in real life

totalN = len(result)
print(totalN)

def downloader(counter=0, default_decounter=250):
    """Function used to download POSCAR files and fill up datalist with information from Aflow
    Works with search result, has two parameters:
    counter - The search entry where we start downloading from (technically we always want zero, only added for flexebility)
    decounter - number of entries we download in a batch before waiting for some time in order not to overburden aflow with requests"""

    with open('./datalist.csv', 'a') as f:        # we create a csv file to write info into
        f.write("ID,compound,lattice_system,spacegroup,volume_cell,moment_cell,mag_field,mag_sites,comment1,comment2")
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