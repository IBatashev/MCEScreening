import time
import random
import pandas as pd
from shutil import move
import os
import subprocess

datadir = './datadir/'
workdir = './workdir'

# to run = files in dir
# complete = no dir
# fail = fail file in ID dir
# in progress = no files inside dir
# skip = no dir

# we have a inputdir where all prepared input files are and wich we modify as noted above,
# a workdir to where files are moved during calculaion
# a outdir where results of calculations are moved to after archivation is done

# a separate tool function that will make up a progress report upon request (uses datalist as starting
# list and fills up status table by parsing inputdir with status of each entry determined as listed above).

################################################################################################

def do(job_type):
    source = datadir + str(item) + '/' + job_type
    files = os.listdir(source)
    if files:
        for f in files:
            move(source + '/' + f, workdir)
        df.loc[item, job_type] = '1111'
        df.to_csv('status.csv')

        callvasp = ''
        ########
        # p1 = subprocess.Popen([callvasp], stdout=subprocess.PIPE)
        # output = p1.communicate()[0]
        ########
        output = 'blabla\n bla reached required accuracy bla' # FOR TESTS REMOVE
        ##########
        with open('./vasprun.out', 'a') as f:  # we create an output file with vasp response
            f.write(output)
        success_line = 'reached required accuracy'
        if success_line in output:
            df.loc[item, job_type] = 'DONE'
            # archive
            # move to savedir
            # clean
            return "complete"
        else:
            df.loc[item, job_type] = 'FAIL'
            # archive
            # move to savedir
            # clean
            return "failed"

    if not files:
        return "no files"


time.sleep(random.randrange(1,10)/2)                                  # sleep for random time 1 - 5 seconds with 0.5 second difference

df = pd.read_csv('status.csv', index_col=0, sep='\t', dtype = str)    # open the status file containing our database entries
result = 'no tasks remain'

for item in df.index.tolist():
    if df.loc[item, 'undeformed'] == '0000':
        result = str(item) + ' undeformed ' +  do('undeformed')
        break
    elif df.loc[item, 'undeformed'] == 'DONE':
        if df.loc[item, 'a'] == '0000':
            result = str(item) + ' a ' +  do('a')
            break
        elif df.loc[item, 'b'] == '0000':
            result = str(item) + ' b ' + do('b')
            break
        elif df.loc[item, 'c'] == '0000':
            result = str(item) + ' c ' + do('c')
            break
        else:
            pass

print(result)
