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

# ------------------------------------------------------------------------------------------------------- #
#  _____                                           _       _____ _             _     _   _                #
# /  __ \                                         | |     /  ___| |           | |   | | | |               #
# | /  \/ ___  _ __ ___  _ __ ___   __ _ _ __   __| |___  \ `--.| |_ __ _ _ __| |_  | |_| | ___ _ __ ___  #
# | |    / _ \| '_ ` _ \| '_ ` _ \ / _` | '_ \ / _` / __|  `--. \ __/ _` | '__| __| |  _  |/ _ \ '__/ _ \ #
# | \__/\ (_) | | | | | | | | | | | (_| | | | | (_| \__ \ /\__/ / || (_| | |  | |_  | | | |  __/ | |  __/ #
#  \____/\___/|_| |_| |_|_| |_| |_|\__,_|_| |_|\__,_|___/ \____/ \__\__,_|_|   \__| \_| |_/\___|_|  \___| #
#                                                                                                         #
# ------------------------------------------------------------------------------------------------------- #


trashdir = 'X:/MCES/MP/trash'
wdatalist = 'X:/MCES/MP/fail_cn07_introduction_error.csv'
outdir = 'X:/MCES/MP/outdir'

outdir_remove(wdatalist, outdir, trashdir)
