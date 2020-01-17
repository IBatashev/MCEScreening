# MCEScreening
Screening for high-performance magnetocaloric materials using information from crystallography databases.

Written in Python3

Essentially MCEScreening consists of three parts:
1. _*Screener*_ - used to download data from Aflow and apply various screening criteria.   
2. _*Creator*_ - a set of python scripts to create VASP input files for all calculations. Can be run locally (and transfer resulting file structure to cluster) or on cluster directly.
3. _*Executor*_ - A bash script submited to queueing system to execute each calculation.

Notes
---
- Moments for MAGMOM are now chosen from Slater curves for magnetic elements (3d and 4d) all other elements set to 1.0 (maybe change to 0.5 - need to check after test runs)
- Potentials from PBE .54 (according to recommendations from vaspwiki). ENCUT set as ENMAX*1.3
- Magnetic field is calculated for PRIMITIVE lattice not the conventional one used in calculations, still the result should hold

Database ChangeLog
---
| # | Date Created | Number of entries | Comment |
|---|---|---|---|
|1|23.09.19|8044|As of now mostly useless. Downloaded in 3 parts and then manually combined, archive also contains initial python scripts and some ‘in progress’ files. |
|2|04.12.19|8970|Relatively adequate second attempt. Reworked list of allowed elements – but still with radioactive...Unfortunately, missing python script used for downloading it.|
|3|11.12.19| 28 |Small test subset made from 2 with python script included in archive. Random seed  = 1|
|4|22.12.19|8970|More files now – separate folders including aflow files (which contain all possible non-file tags). Instead of CONTCARs we now have aflow_structure files for both before and after relaxation containg a LOT of structural info. Contains python script used for downloading. Includes python script used for downloading. Aflow files were downloaded a week |
|5|07.01.19|8603|Updated version of 4, radioactive things now removed. Includes script that was used to remove entries containing  U, Po, Th|
TODO
---
- Ask Gilles about W_pv and At_d in PBE .54
- Redo POSCAR_maker (and all other that were based on POSCAR file from aflow) - write lattice matrix from lattice parameters and symmetry instead of copying from original POSCAR
- Add check to make sure direct coordinates are supplied in initial POSCAR from aflow (POSCAR_maker)
- ~~add ase package to project directory~~ with new info from database I no longer really need ase 
- Fix lattice recognition for making deformations
- add error messages/exceptions everywhere so we can troubleshoot
- Add scripts for auotomatic processing - see 'custodian' package
---
    nn and nnn distance as creening parameter - check
    shifting magnetic sublattices against each other?
