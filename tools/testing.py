from pymatgen.core import Structure, IStructure
from pymatgen.symmetry.analyzer import SpacegroupAnalyzer
from pymatgen.ext.matproj import MPRester

from pymatgen.io.vasp.inputs import Poscar

##mpid = 'mp-148'
##
##rester = MPRester()
##structure = rester.get_structure_by_material_id(mpid)
##print(structure)

structure = IStructure.from_file('POSCAR')
print(structure)

a = SpacegroupAnalyzer(structure, symprec=1e-4)
structure = a.get_conventional_standard_structure()
# print(structure)
print(structure.get_space_group_info())

p = Poscar(structure)
print(p)

p.write_file('POSCAR.pym.conv')

structure = a.get_primitive_standard_structure()
# print(structure)
print(structure.get_space_group_info())


p = Poscar(structure)
print(p)

p.write_file('POSCAR.pym.prim')

