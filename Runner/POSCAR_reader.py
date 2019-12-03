def read(poscar_file):
    """Takes POSCAR file as input and returns it's content as list of strings"""

    poscar = open(poscar_file, "r")
    poscar_content = []
    for line in poscar:
        poscar_content.append(line)
    poscar.close()
    return poscar_content