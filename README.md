# MCEScreening
Screening for high-performance magnetocaloric materials using information from crystallography databases.

Written in Python3

Essentially MCEScreening consists of three parts:
1. _*Downloader*_ - used to download data from various material databases.
2. _*Screener*_   a set of python scripts to apply various screening criteria before and after DFT calculations.   
3. _*Creator*_ - a set of python scripts to create VASP input files for all calculations. Can be run locally (and transfer resulting file structure to cluster) or directly on the cluster .
4. _*Executor*_ - A script submited to queueing system to execute a single calculation. 

File Structure
---

1. datalist - a listing of all entries with some of their properties. Created during initial download of the database. We upate it with new parameters and use it for screening.
2. datadir - contains files from Aflow as downloaded. There are now 4 files for each entry: INCAR, structure before relaxation, structure after relaxation and a file with all entry properties that are provided by aflow. (only relaxed structure and aflow file are really necessary/important).
3. inputdir - contains files required for performing all the calulations. This dir is prepared in advance by creator.py and then moved to supercomputer. It is removed emptied during the calculations.
4. outdir - cointains results of calculations for all entries (zipped). This dir is created during the calculations.
5. statuslist - record of progress of computations. Prepared together with inputdir and updated on request during or after calculation based on content of inputdir and outputdir.


