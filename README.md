# MCEScreening
Screening for high-performance magnetocaloric materials using information from crystallography databases.

TODO:
- Test if sym_detector works as intended, if not lower symmetry tolerance?
- Add scripts for auotomatic processing - see 'custodian' package
- Fix the thing about volume increase constant in initial POSCAR from aflow
- Add check to make sure direct coordinates are supplied in initial POSCAR from aflow
- All latice deformations have to be properly tested
- add error messages/exceptions everywhere so we can troubleshoot
- download 8000? files again and sort them. (aflow_ICSD)
- create a semi-random set of  materials to run quck tests
- fix rhombohedral lattice deformation

Essentially MCEScreening consists of two parts:
1. _Screener_ - used to download data from Aflow and apply various screening criteria.   
2. _Runner_ - a set of scripts used to run various vasp calculations. These work on cluster nodes, and preferably only 
have a very basic set of dependencies to avoid the need to insatal many extra packages.

