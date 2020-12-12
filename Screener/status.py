import os
import pandas as pd
import tqdm


def status_before(datalist, inputdir, start=-1):
    """ Function that reads through inputdir and lists undeformed structure and all it's
    deformations and writes result to a new .csv file."""
    counter = 0
    df = pd.read_csv(datalist, index_col=0, sep=',')

    with tqdm.tqdm(total=len(df.index)) as pbar:        # A wrapper that creates nice progress bar
        pbar.set_description("Processing datalist")
        for item in df.index.tolist():
            pbar.update(1)                              # Updating progress bar at each step
            path = inputdir + '/' + str(item)
            deformations = os.listdir(path)
            for deformation in deformations:
                df.loc[item, str(deformation)] = 'to_run'
                counter = counter+1
                df.to_csv((datalist.replace('.csv', '_beforeRun.csv')))
    print('Total number of calculations to perform ', counter)


def status_after(datalist, outdir):
    """A function to display progress of screening or prepare a report on completed screening - a csv file with all IDs
    their with their deformations and number of warnings failures"""

    sub_calculations_list = ['a_inc', 'a_dec', 'b_inc', 'b_dec', 'c_inc', 'c_dec']
    fail_counter = 0
    run_counter = 0
    total_warning_counter = 0
    df = pd.read_csv(datalist, index_col=0, sep=',')

    with tqdm.tqdm(total=len(df.index)) as pbar:        # A wrapper that creates nice progress bar
        pbar.set_description("Processing datalist")
        for item in df.index.tolist():
            pbar.update(1)                              # Updating progress bar at each step

            warn_counter = 0
            path = outdir + '/' + str(item)
            if os.path.exists(path):                                    # first we look if this entry was calculated, if folder in outdir does not exist it means it was not run at all and status will remain 'to_run'

                if df.loc[item, 'undeformed'] == 'to_run':              # Next, we check if undeformed was calculated (this is purely for the case when we use this script after we calculated some specific deformation in a new run and want to merge results)
                    if os.path.exists(path + '/' + 'undeformed'):
                        files = os.listdir(path + '/' + 'undeformed')
                        run_counter = run_counter + 1
                        if 'warning' in files:
                            with open(path + '/' + 'undeformed' + '/warning', 'r') as f:
                                for line in f:
                                    warn_counter = warn_counter + 1
                                    total_warning_counter = total_warning_counter + 1
                        if 'fail' in files:                                 # if there is a fail flag in undeformed subfolder all deformations that were supposed to run for this compounds are marked as failures
                            df.loc[item, 'undeformed'] = 'failed'
                            with open(path + '/' + 'undeformed' + '/fail', 'r') as failfile:
                                df.loc[item, 'undeformed_fail_reason'] = failfile.readline().strip('\n')
                            fail_counter = fail_counter + 1
                            for sub_calc in sub_calculations_list:
                                if df.loc[item, sub_calc] == 'to_run':
                                    df.loc[item, sub_calc] = 'failed'
                                    fail_counter = fail_counter + 1
                            continue                                        # we no longer need to separately check subruns anymore so we skip to next entry
                        else:
                            df.loc[item, 'undeformed'] = 'completed'
                    else:
                        df.loc[item, 'undeformed'] = 'failed'
                        df.loc[item, 'undeformed' + '_fail_reason'] = 'not calculated'
                        fail_counter = fail_counter + 1
                        for sub_calc in sub_calculations_list:
                            if df.loc[item, sub_calc] == 'to_run':
                                df.loc[item, sub_calc] = 'failed'
                                fail_counter = fail_counter + 1

                for sub_calc in sub_calculations_list:
                    if df.loc[item, sub_calc] == 'to_run':
                        if os.path.exists(path + '/' + sub_calc):
                            files = os.listdir(path + '/' + sub_calc)
                            run_counter = run_counter + 1
                            if 'warning' in files:
                                with open(path + '/' + sub_calc + '/warning', 'r') as f:
                                    for line in f:
                                        warn_counter = warn_counter + 1
                                        total_warning_counter = total_warning_counter + 1
                            if 'fail' in files:
                                df.loc[item, sub_calc] = 'failed'
                                with open(path + '/' + sub_calc + '/fail', 'r') as failfile:
                                    df.loc[item, sub_calc+'_fail_reason'] = failfile.readline().strip('\n')
                                fail_counter = fail_counter + 1
                            else:
                                df.loc[item, sub_calc] = 'completed'
                        else:
                            df.loc[item, sub_calc] = 'failed'
                            df.loc[item, sub_calc + '_fail_reason'] = 'not calculated'
                            fail_counter = fail_counter + 1
                df.loc[item, 'warnings'] = warn_counter

    df.to_csv((datalist.replace('.csv', '_afterRun.csv')))
    print('Total number of calculations', run_counter)
    print('Total number of failures', fail_counter)
    print('Total number of warnings', total_warning_counter)


def make_inputdir_for_rerun(datalist, inputdir_initial, inputdir_rerun):
    """Takes datalist, looks what runs failed and takes corresponding folders from initial inputdir
    to create an inputdir to use for rerun"""


def separate_failed_entries():
    """Moves all entries in failed.csv to a separate folder from outdir (and sorts them by fail type?)"""


def sieve_for_success(datalist):
    """Separates datalist into one for fully sucessful entries and one for entries with some kind of failures
    We actually only need it if we want to use Screener module on a run that was not completely successful yet, as
    Screener cannot properly deal with entries that have failures. .csv of failures is just an extra for convenience"""

    success_counter = 0
    fail_counter = 0
    df = pd.read_csv(datalist, index_col=0, sep=',')
    df_fail = df
    for item in df.index.tolist():
        if (df.loc[item, 'undeformed'] != 'failed') and \
                (df.loc[item, 'a_dec'] != 'failed') and \
                (df.loc[item, 'a_inc'] != 'failed') and \
                (df.loc[item, 'b_dec'] != 'failed') and \
                (df.loc[item, 'b_inc'] != 'failed') and \
                (df.loc[item, 'c_dec'] != 'failed') and \
                (df.loc[item, 'c_inc'] != 'failed'):
            success_counter = success_counter + 1
            df_fail = df_fail.drop([item], axis=0)
        else:
            df = df.drop([item], axis=0)
            fail_counter = fail_counter + 1

    print('Number of successful entries', success_counter)
    print('Number of failed entries', fail_counter)

    df.to_csv(datalist.replace(".csv", '_success_sieved.csv'))
    df_fail.to_csv(datalist.replace(".csv", '_failed_sieved.csv'))


def sieve_for_torun(datalist):

    success_counter = 0
    torun_counter = 0
    df = pd.read_csv(datalist, index_col=0, sep=',')
    df_torun = df
    for item in df.index.tolist():
        if (df.loc[item, 'undeformed'] != 'to_run') and \
                (df.loc[item, 'a_dec'] != 'to_run') and \
                (df.loc[item, 'a_inc'] != 'to_run') and \
                (df.loc[item, 'b_dec'] != 'to_run') and \
                (df.loc[item, 'b_inc'] != 'to_run') and \
                (df.loc[item, 'c_dec'] != 'to_run') and \
                (df.loc[item, 'c_inc'] != 'to_run'):
            success_counter = success_counter + 1
            df_torun = df_torun.drop([item], axis=0)
        else:
            df = df.drop([item], axis=0)
            torun_counter = torun_counter + 1

    print('Number of complete entries', success_counter)
    print('Number of not yet run entries', torun_counter)

    df.to_csv(datalist.replace(".csv", '_complete.csv'))
    df_torun.to_csv(datalist.replace(".csv", '_torun_sieved.csv'))

# ------------------------------------------------------------------------------------------------------- #
#  _____                                           _       _____ _             _     _   _                #
# /  __ \                                         | |     /  ___| |           | |   | | | |               #
# | /  \/ ___  _ __ ___  _ __ ___   __ _ _ __   __| |___  \ `--.| |_ __ _ _ __| |_  | |_| | ___ _ __ ___  #
# | |    / _ \| '_ ` _ \| '_ ` _ \ / _` | '_ \ / _` / __|  `--. \ __/ _` | '__| __| |  _  |/ _ \ '__/ _ \ #
# | \__/\ (_) | | | | | | | | | | | (_| | | | | (_| \__ \ /\__/ / || (_| | |  | |_  | | | |  __/ | |  __/ #
#  \____/\___/|_| |_| |_|_| |_| |_|\__,_|_| |_|\__,_|___/ \____/ \__\__,_|_|   \__| \_| |_/\___|_|  \___| #
#                                                                                                         #
# ------------------------------------------------------------------------------------------------------- #


# datalist_before = 'D:/MCES/aflow/errors_run1/datalist.csv'
# datalist_after = 'D:/MCES/aflow/errors_run1/datalist_beforeRun.csv'
# outdir = 'D:/MCES/aflow/errors_run1/outdir'
# inputdir = 'D:/MCES/aflow/errors_run1/inputdir'

# outdir = 'D:/MCES/MP/batch1/outdir'

inputdir = 'D:/MCES/MP/inputdir'
datalist_before = ''

# datalist_before = 'D:/MCES/MP/datalist_lattfix_updated_sieved.mag.field_sieved.mag.sites_no.duplicates.csv'

outdir = 'D:/MCES/MP/outdir'
datalist_after = 'D:/MCES/MP/datalist_lattfix_updated_sieved.mag.field_sieved.mag.sites_no.duplicates_beforeRun.csv'


# status_before(datalist_before, inputdir)

# status_after(datalist_after, outdir)


# sieve_for_success('D:/MCES/MP/datalist_lattfix_updated_sieved.mag.field_sieved.mag.sites_no.duplicates_beforeRun_afterRun.csv')
# sieve_for_torun('D:/MCES/MP/datalist_lattfix_updated_sieved.mag.field_sieved.mag.sites_no.duplicates_beforeRun_afterRun_success_sieved.csv')
# status_before('MnAs.csv', 'MnAs_in')
status_after('MnAs_beforeRun.csv', 'MnAs_o')