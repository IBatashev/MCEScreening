# to combine several runs of the same database in one filestructure e.g. add BEXT run results to previously done magnetoelastic run
import os
import pandas as pd
import shutil


def outdir_remove(datalist, outdir, trash=None):
    """Removes entries from outdir that are listed in supplied datalist.
    Optionally provide trash directory path, then files would be moved there instead of being straight up deleted"""

    df = pd.read_csv(datalist, index_col=0, sep=',')
    for item in df.index.tolist():
        if trash is None:
            os.remove(outdir + '/' + str(item))
        else:
            shutil.move(outdir + '/' + str(item), trash + '/' + str(item))


def outdir_copy(datalist, source, destination):
    """"Copies output folders with ids matching ids in the datalist to a separate directory selected by user"""
    df = pd.read_csv(datalist, index_col=0, sep=',')
    for item in df.index.tolist():
        shutil.copytree(source + '/' + str(item), destination + '/' + str(item))


def file_copy(datalist, source, destination):
    """"Copies files with ids matching ids in the datalist to a separate directory selected by user"""

    df = pd.read_csv(datalist, index_col=0, sep=',')
    for item in df.index.tolist():
        dest_fpath = destination + '/' + str(item) + '/' + 'sym'
        os.makedirs(dest_fpath, exist_ok=True)
        shutil.copyfile(source + '/' + str(item) + str('.cif'), dest_fpath + '/' + str(item) + str('.cif'))


def rename_subdir(datalist, outdir, before, after):
    """Renames specified subfolders for entries from dir that are listed in supplied datalist.
    Optionally provide trash directory path, then files would be moved there instead of being straight up deleted"""
    df = pd.read_csv(datalist, index_col=0, sep=',')
    for item in df.index.tolist():
        try:
            shutil.move(outdir + '/' + str(item) + '/' + str(before), outdir + '/' + str(item) + '/' + str(after))
        except:
            print(item)


def outdir_remove_subdir(datalist, outdir, subdir, trash=None, ):
    """Removes specified subfolders for entries from dir that are listed in supplied datalist.
    Optionally provide trash directory path, then files would be moved there instead of being straight up deleted"""

    df = pd.read_csv(datalist, index_col=0, sep=',')
    for item in df.index.tolist():
        try:
            if trash is None:
                os.remove(outdir + '/' + str(item) + '/' + str(subdir))
            else:
                shutil.move(outdir + '/' + str(item) + '/' + str(subdir), trash + '/' + str(item) + '/' + str(subdir))
        except:
            print('some problem with', item)
# ------------------------------------------------------------------------------------------------------- #
#  _____                                           _       _____ _             _     _   _                #
# /  __ \                                         | |     /  ___| |           | |   | | | |               #
# | /  \/ ___  _ __ ___  _ __ ___   __ _ _ __   __| |___  \ `--.| |_ __ _ _ __| |_  | |_| | ___ _ __ ___  #
# | |    / _ \| '_ ` _ \| '_ ` _ \ / _` | '_ \ / _` / __|  `--. \ __/ _` | '__| __| |  _  |/ _ \ '__/ _ \ #
# | \__/\ (_) | | | | | | | | | | | (_| | | | | (_| \__ \ /\__/ / || (_| | |  | |_  | | | |  __/ | |  __/ #
#  \____/\___/|_| |_| |_|_| |_| |_|\__,_|_| |_|\__,_|___/ \____/ \__\__,_|_|   \__| \_| |_/\___|_|  \___| #
#                                                                                                         #
# ------------------------------------------------------------------------------------------------------- #



wdatalist = 'D:/MCES/COD/step3trimclort.csv'
trashdir = 'D:/MCES/COD/trash'
outdir = 'D:/MCES/COD/outdir'
# outdir2 = 'D:/MCES/MP/outdir_V/done'


# outdir = 'D:/MCES/COD/datadir'



# rename_subdir(wdatalist, outdir, 'a_1.16', 'a_1.05')
# outdir_remove_subdir(wdatalist, outdir, 'c_0.95', trashdir)
#outdir_copy('D:/MCES/MP/more_points_without_O.csv', 'D:/MCES/MP/outdir', 'D:/MCES/MP/new_outdir')
file_copy('D:/MCES/MP/datalist_for_Tc_less_70_atoms3.csv', 'D:/MCES/MP/datadir', 'D:/MCES/MP/inputdir3_Tc')
