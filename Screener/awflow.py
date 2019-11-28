from aflow import *
import pandas
import numpy as np
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


# df = pandas.read_csv('data.csv')
# df.to_csv('hrdata_modified.csv')

totalN = len(result)
print(totalN)
ds = np.empty([0,8])

###############################################################################################################
# GRABBER THAT REALLY WORKS
# counter = 1599
# counter = 2096
# counter = 2563
# counter = 3064
# counter = 3218
# counter = 3677
# counter = 4108
# counter = 7079
# counter = 7360
# decounter = 250

while counter <= totalN:
    result[counter].files["CONTCAR.relax.vasp"]("./data/"+str(counter))
    newrow = str(counter)+','+ str(result[counter].compound)+','+str(result[counter].lattice_system_relax.strip('\n'))+','+ str(result[counter].spacegroup_relax)+','+str(result[counter].volume_cell)+','+ str(result[counter].spin_cell)+','+'0,0\n'
    with open('./data.txt', 'a') as f:
        f.write(newrow)
    if decounter == 0:
        print(counter, '   ', (counter/totalN)*100, '% done')
        time.sleep(300)
        decounter = 250
    counter = counter + 1
    decounter = decounter - 1

#####################################################################################################################################################################################
# modified version I used to get some missing data
# counter = 2914
# decounter = 500
# while counter <= 3064:
#     newrow = str(counter)+','+ str(result[counter].compound)+','+str(result[counter].lattice_system_relax.strip('\n'))+','+ str(result[counter].spacegroup_relax)+','+str(result[counter].volume_cell)+','+ str(result[counter].spin_cell)+','+'0,0\n'
#     with open('./data2.txt', 'a') as f:
#         f.write(newrow)
#     if decounter == 0:
#         print(counter, '   ', (counter/totalN)*100, '% done')
#         time.sleep(300)
#         decounter = 500
#     counter = counter + 1
#     decounter = decounter - 1

