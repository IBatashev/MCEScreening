# Step 1. Downloader

Download structure files from the database of choice using _Downloader_ scripts (Awflow_downloader.py, COD_downloader.py, Materials_Project_downloader.py and ICSD_get_ids.py + ICSD_get_cifs_by_ids.py)
. Initial screening parameters are hard-coded in the scripts - modify code to change them.

<ins>Input</ins>: None

<ins>Resulting files</ins>: *datalist.csv* and *datadir* folder with .cif files for all structures

Notes:
   - Aflow is likely outdated
   - Material Project will soon migrate to new API and will be outdated (check database website)
   - COD has to be fully downloaded from https://wiki.crystallography.net/howtoobtaincod/
   - ICSD scripts need to run on cn52 or any machine with access to ICSD SQL database
   - Alternatively simply create  *datalist.csv* and *datadir* manually (for tests or smaller scale runs)


# Step 2. Screener 
Apply screening criteria that are available just from structure.
This is database specific step use either *screener.py* (Aflow, MP), *screener_COD.py* or *screener_ICSD.py*

<ins>Input</ins>: None

<ins>Resulting files</ins>: *datalist.csv* and *datadir* folder with .cif files for all structures


# Step 3. Creator
Create input files for DFT calculations based on selected  *datalist.csv* and *datadir*. 
Main script *creator.py* contains various helper functions and will also call for *INCAR_maker.py*, *POSCAR_maker.py*, *POTCAR_maker.py* scripts.
Make choice of what to do at the bottom of *creator.py* script after #COMMANDS START HERE# comment line. Here you need to
select paths for *datalist.csv*,  *datadir*, *inputdir*(where VASP inputs will be created) and call function 

creator(wdatalist, wdatadir, def_factor=0.05, hydrostatic=True, uniaxial=True, undeformed=True, field=False)

Notes:
- def_factor selects deformation % for both hydrostatic and uniaxial types of deformations
- undeformed - no adjustments to lattice or field (normal VASP calculation)
- hydrostatic=True will modify whole volume 
- uniaxial will deform axis individually
- field activates external field in VASP (BEXT)



# Step 3.5. Screener
Index input files created in step 2 to track progress of calculations. This is optional but highly recommended.