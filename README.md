# MCEScreening
Screening for high-performance magnetocaloric materials using information from crystallography databases.

Written in Python3

Essentially MCEScreening consists of three parts:
1. _*Screener*_ - used to download data from Aflow and apply various screening criteria.   
2. _*Creator*_ - a set of python scripts to create VASP input files for all calculations. Can be run locally (and transfer resulting file structure to cluster) or directly on the cluster .
3. _*Executor*_ - A script submited to queueing system to execute a single calculation. 

File Structure
---

1. datalist - a listing of all entries with some of their properties. Created during initial download of the database. We upate it with new parameters and use it for screening.
2. datadir - contains files from Aflow as downloaded. There are now 4 files for each entry: INCAR, structure before relaxation, structure after relaxation and a file with all entry properties that are provided by aflow. (only relaxed structure and aflow file are really necessary/important).
3. inputdir - contains files required for performing all the calulations. This dir is prepared in advance by creator.py and then moved to supercomputer. It is removed emptied during the calculations.
4. outdir - cointains results of calculations for all entries (zipped). This dir is created during the calculations.
5. statuslist - record of progress of computations. Prepared together with inputdir and updated on request during or after calculation based on content of inputdir and outputdir.


Database ChangeLog
---
| # | Date Created | Number of entries | Comment |
|---|---|---|---|
|a1|23.09.19|8044|As of now mostly useless. Downloaded in 3 parts and then manually combined, archive also contains initial python scripts and some ‘in progress’ files. |
|a2|04.12.19|8970|Relatively adequate second attempt. Reworked list of allowed elements – but still with radioactive...Unfortunately, missing python script used for downloading it.|
|a3|11.12.19| 28 |Small test subset made from 2 with python script included in archive. Random seed  = 1|
|a4|22.12.19|8970|More files now – separate folders including aflow files (which contain all possible non-file tags). Instead of CONTCARs we now have aflow_structure files for both before and after relaxation containg a LOT of structural info. Contains python script used for downloading. Includes python script used for downloading. Aflow files were downloaded a week |
|a5|07.01.19|8603|Updated version of 4, radioactive things now removed. Includes script that was used to remove entries containing  U, Po, Th|
|a6|09.02.20| 31 |Small test subset made from 5 using make_test_db.py with random seed = 1. After random generation manually added Fe2P(ID=6295), FeRh(ID=6373) and LaFeSi(ID=5565)|
|a7|31.03.20| 20 |Small test subset made from 5 using make_test_db.py with random seed = 1. All compounds are hexagonal (Fe2P(ID=6295) was already in from random|
|MP1|13.05.20|12160| First snapshot from Materials Project same general settings as with aflow but, more info and structure files are in .cif format|

Run Results
---
| # | Date Performed | Comment | Result Comment
|---|---|---|---|
|1|11.02.20|TestDB (database 6) with deformation constant of 1.2 (20% increasae)| Speed ok; no errors; some entries lost moments - need to check; need to do different deformation constant (0.8) maybe several steps

Notes
---
- Moments for MAGMOM are now chosen from Slater curves for magnetic elements (3d and 4d) all other elements set to 1.0 (maybe change to 0.5 - need to check after test runs)
- Potentials from PBE .54 (according to recommendations from vaspwiki, except for Mg, W_pv and At_d which have more recent updates). ENCUT set as ENMAX*1.3

TODO
---
- (optional) Add scripts for automatic processing of the calculations - see 'custodian' package
- clean all project files from test things, make them easier to use - less variables to change before runs
- add instructions on how to use MCES add list of steps and in what order to do them and also write which files are created at each step

- new flowchart for executor
- Gd5Ge2Si2 MCL reference
- add rounding to numbers before they go into .csv
- wtf is wrong with entries 2870 and 7513 why are the structures not recognized by pymatgen?
- Add function to screener_after to parse messages from warnlist and badlist and add them to datalist somehow - just write number of warnings? and failures?
- add check to screener_after to  warn about big change in mag field after calculation compared to before
- make a script for screen after that gathers result from the run - failures and reasons etc to make a csv file with successfully completed compounds that we can further analyse
- update readme
- make option for creator to include bext
- write folder joiner tool
- rewrite function description to display parameter tips
- simplify executor

- Collector for results:
    - symmetries (bravais + group number) + lattice parameters from OUTCAR - ALAT and C/A, B/A
    - moments (total)
    - moments by elements
    - geometry (volume + a,b,c) 
    - number of atoms in structure
    - if the deformations are done the way we expected
    - average memory use
    - total number of warnings

---
    nn and nnn distance as screening parameter - check
    how number (distances) of nn changes in FM and AFM state - not nece, requires relaxatino   
    free volume as screening parameter (compare total cell volume with volume occupied by atoms (covalent radius))?)
    shifting magnetic sublattices against each other?
---

1st make new inputdir with creator
2nd create a .csv from this inputdir to know which deformation we are going to performance
3rd run VASP calculations
4th check which runs failed and create a new .csv that only contains failures
5th run creator to make new inputdir from failures.csv
6th run VASP calculations now we want to actually retry failed runs with adjusted incars

7th join run results in outdir, replacing failed subfolders and adding subfolders that were not previously there
