# MCEScreening
Screening for high-performance magnetocaloric materials using information from crystallography databases.

Written in Python3

Essentially MCEScreening consists of two parts:
1. _*Screener*_ - used to download data from Aflow and apply various screening criteria.   
2. _*Runner*_ - a set of python scripts to create VASP input files for all calculations. These work on cluster nodes, 
and preferably only have a very basic set of dependencies to avoid the need to insatal many extra packages.
A bash script is then submited to queueing system to execute each calculation.

Notes
---
- Moments for MAGMOM are now chosen from Slater curves for magnetic elements (3d and 4d) all other elements set to 1.0 (maybe change to 0.5 - need to check after test runs)
- Potentials from PBE .54 (according to recommendations from vaspwiki). ENCUT set as ENMAX*1.3

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
- Perform initial screening - calculate internal fileds and sites
- Ask Gilles about W_pv and At_d in PBE .54
- Figure out what elements should be in mag_site calculator magnetic list
- Problem with structure files - ABC... instead of elements, mag sites is affected and other POSCAR-related files
- Redo POSCAR_maker (and all other that were based on POSCAR file from aflow) - write lattice matrix from lattice parameters and symmetry instead of copying from original POSCAR
- Add check to make sure direct coordinates are supplied in initial POSCAR from aflow (POSCAR_maker)
- Test if mag_sites works as intended, if not lower symmetry tolerance? Also no longer need to employ ase or phonopy - just grab wyckof from edata.relax.out
- ~~add ase package to project directory~~ with new info from database I no longer really need ase 
- Fix lattice recognition
- Get doi from structure file on line starting with "NO RECURSION ... doi:...""
- add error messages/exceptions everywhere so we can troubleshoot
- Add scripts for auotomatic processing - see 'custodian' package
---
    nn and nnn distance as creening parameter - check
    shifting magnetic sublattices against each other?
