import os
import pandas as pd
import tqdm



def status_before(datalist, calc_dict):
    """ Function that reads through datalists and determines what calculations have to be run,
    depending on symmetry and supplied list of calculation types
     writes a new .csv file."""

    counter = 0
    df = pd.read_csv(datalist, index_col=0, sep=',')

    if 'undeformed' in list(calc_dict):
        df['undeformed'] = 'to_run'
        counter = counter + len(df.index)
    if 'Applied_Field' in list(calc_dict):
        df['Applied_Field'] = 'to_run'
        counter = counter + len(df.index)
    if 'volumetric' in list(calc_dict):
        for i in calc_dict['volumetric']:
            df['V_'+str(i)] = 'to_run'
            counter = counter + len(df.index)

    with tqdm.tqdm(total=len(df.index)) as pbar:        # A wrapper that creates nice progress bar
        pbar.set_description("Processing datalist")    # the speed can be vastly improved here. df.loc is NOT a good way to do this operation, but no time to fix now
        for item in df.index.tolist():
            pbar.update(1)                              # Updating progress bar at each step
            symmetry = df.loc[item, 'lattice_system']
            if symmetry == 'cubic' or symmetry == 'rhombohedral' or symmetry == 'trigonal':
                par_list = ['a_']
            if symmetry == 'hexagonal' or symmetry == 'tetragonal':
                par_list = ['a_', 'c_']
            if symmetry == 'monoclinic' or symmetry == 'orthorhombic' or symmetry == 'triclinic':
                par_list = ['a_', 'b_', 'c_']
            if 'uniaxial' in list(calc_dict):
                for i in calc_dict['uniaxial']:
                    calc_list = [x+i for x in par_list]
                    for j in calc_list:
                        df.loc[item, str(j)] = 'to_run'
                        counter = counter + 1

    df.to_csv((datalist.replace('.csv', '_beforeRunning.csv')))
    print('Total number of calculations to perform ', counter)


def calc_list(calc_dict):
    sub_calculations_list = []
    if 'undeformed' in list(calc_dict):
        sub_calculations_list.append('undeformed')
    if 'Applied_Field' in list(calc_dict):
        sub_calculations_list.append('Applied_Field')
    if 'volumetric' in list(calc_dict):
        for i in calc_dict['volumetric']:
            sub_calculations_list.append('V_' + str(i))
    if 'uniaxial' in list(calc_dict):
        for i in calc_dict['uniaxial']:
            par_list = ['a_', 'b_', 'c_']  # we expect to always have a set of calculations with all symmetry types, even if not, does not matter
            calc_list = [x + i for x in par_list]
            for j in calc_list:
                sub_calculations_list.append(str(j))
    return sub_calculations_list


def status_after(datalist, outdir, calc_dict):
    """A function to display progress of screening or prepare a report on completed screening - a csv file with all IDs
    their with their deformations and number of warnings failures"""
    sub_calculations_list = calc_list(calc_dict)
    critical_failure_count = 0
    fail_counter = 0
    run_counter = 0
    total_warning_counter = 0
    df = pd.read_csv(datalist, index_col=0, sep=',')

#    critical_failure if either undeformed or all paired deformations not performed

    with tqdm.tqdm(total=len(df.index)) as pbar:        # A wrapper that creates nice progress bar
        pbar.set_description("Processing datalist")
        for item in df.index.tolist():
            pbar.update(1)                              # Updating progress bar at each step

            critical_failure = False

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
                            critical_failure_count = critical_failure_count + 1
                            critical_failure = True
                            df.at[item, 'crit_fail'] = 'Yes'
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
            else:
                critical_failure_count = critical_failure_count + 1
                critical_failure = True

            # now we decide which rows can be further processed
            df2 = df.loc[item, sub_calculations_list]
            df2.dropna(how='all', inplace=True)

            df2 = df2.to_frame()
            df2['type'] = df2.index
            df2.columns = ['result', 'type']
            df2['type'] = df2['type'].str[2:]

            if df2.loc['undeformed', 'result'] == 'completed':
                latlist = list(df2.index)
                try:
                    latlist.remove('undeformed')
                    latlist.remove('Applied_Field')
                except: pass

                latlist = set([w[:1] for w in latlist])
                def_failures = False
                for lat in latlist:
                    df3 = df2.loc[df2.index.str.startswith(lat)].copy()
                    df3['type'] = abs(1 - (df3['type']).astype('float')).round(3)
                    df3 = df3.drop(df3[df3.result == 'failed'].index)
                    if len(df3.index) == 0:
                        pass
                    else:
                        pairs = len(df3['type']) - len(df3.drop_duplicates(subset=['type']))
                        if pairs == 0:
                            critical_failure = True
                            def_failures = True
                if def_failures == True:
                    critical_failure_count = critical_failure_count + 1
                    critical_failure = True
            else:
                critical_failure_count = critical_failure_count + 1
                critical_failure = True
            if critical_failure == True:
                df.at[item, 'crit_fail'] = 'Yes'

    df.to_csv((datalist.replace('.csv', '_afterRun.csv')))
    print('Total number of calculations', run_counter)
    print('Total number of failures', fail_counter)
    print('Total number of warnings', total_warning_counter)
    print('Unusable data for ', critical_failure_count, ' compounds')


def sieve_for_success(datalist, calc_dict):
    """Separates datalist into one for fully sucessful entries and one for entries with some kind of failures
    We actually only need it if we want to use Screener module on a run that was not completely successful yet, as
    Screener cannot properly deal with entries that have failures. .csv of failures is just an extra for convenience"""

    columns_to_remove = calc_list(calc_dict)
    fail_reason_columns = [x + '_fail_reason' for x in columns_to_remove]
    columns_to_remove.append(fail_reason_columns)
    columns_to_remove.append('crit_fail')

    df = pd.read_csv(datalist, index_col=0, sep=',')

    # df_suc = df[~df['crit_fail'].isin(['yes']).any(axis=1)]
    try:
        df_suc = df.drop(df[df['crit_fail'] == 'Yes'].index)
        for i in columns_to_remove:
            try:
                df_suc.drop(i, inplace=True, axis=1)
            except:
                pass
    except: df_suc = df

    # df_fail = df[df[calculations_checklist].isin(['failed']).any(axis=1)]
    # df_torun = df[df[calculations_checklist].isin(['to_run']).any(axis=1)]
    df_suc.to_csv(datalist.replace(".csv", '_success_sieved.csv'))
    # df_fail.to_csv(datalist.replace(".csv", '_failed_sieved.csv'))
    # df_torun.to_csv(datalist.replace(".csv", '_torun_sieved.csv'))

    print('Number of successful entries', len(df_suc.index))
    # print('Number of failed entries', len(df_fail.index))
    # print('Number of not yet run entries', len(df_torun.index))


def status_RSPt(datalist, outdir):
    """A function to display progress of RSPt and UppASD screening run and prepare a report - a csv file with current status for all IDs"""
    fail_counter = 0
    success_counter = 0

    run_counter = 0
    df = pd.read_csv(datalist, index_col=0, sep=',')

    with tqdm.tqdm(total=len(df.index)) as pbar:  # A wrapper that creates nice progress bar
        pbar.set_description("Processing datalist")
        for item in df.index.tolist():
            pbar.update(1)  # Updating progress bar at each step

            df.loc[item, 'RSPT_DFT_run'] = 'to_run' # first we assumwe that calculations are not yer complete
            df.loc[item, 'UppASD_run'] = 'to_run'

            path = outdir + '/' + str(item)

            if os.path.exists(path):  # look if this entry was calculated, if folder in outdir does not exist it means it was not run at all and status will remain as 'to_run'
                run_counter = run_counter + 1
                if len(os.listdir(path)) == 0:  # if no files are in the outdir it means that job was cancelled due time limit during Jij calculations or UppASD run
                    df.loc[item, 'RSPT_DFT_run'] = 'failed'
                    df.loc[item, 'RSPT_DFT_run' + '_fail_reason'] = 'cancelled due to time limit after DFT run'
                    df.loc[item, 'UppASD_run'] = 'not calculated'
                    fail_counter = fail_counter + 1
                else:
                    files = os.listdir(path)
                    if 'fail' in files:
                        df.loc[item, 'RSPT_DFT_run'] = 'failed'
                        with open(path + '/' + '/fail', 'r') as failfile:
                            df.loc[item, 'RSPT_DFT_run_fail_reason'] = failfile.readline().strip('\n')
                        fail_counter = fail_counter + 1
                        df.loc[item, 'UppASD_run'] = 'not calculated'
                    else:
                        df.loc[item, 'RSPT_DFT_run'] = 'completed'
                        if 'out_UppASD' in files:  # lastly we check if there is proper output for the UppASD run
                            df.loc[item, 'UppASD_run'] = 'completed'
                            success_counter = success_counter + 1
                            with open(path + '/' + 'out_UppASD', 'r') as UppASD_data:
                                df.loc[item, 'Tc'] = UppASD_data.readlines()[32].split()[5]
                        else:
                            df.loc[item, 'UppASD_run'] = 'not calculated'
                            fail_counter = fail_counter + 1

    df.to_csv((datalist.replace('.csv', '_afterRun.csv')))

    print('Total number of calculations', run_counter)
    print('Number of failed runs', fail_counter)
    print('Number of successful runs', success_counter)


# ------------------------------------------------------------------------------------------------------- #
#  _____                                           _       _____ _             _     _   _                #
# /  __ \                                         | |     /  ___| |           | |   | | | |               #
# | /  \/ ___  _ __ ___  _ __ ___   __ _ _ __   __| |___  \ `--.| |_ __ _ _ __| |_  | |_| | ___ _ __ ___  #
# | |    / _ \| '_ ` _ \| '_ ` _ \ / _` | '_ \ / _` / __|  `--. \ __/ _` | '__| __| |  _  |/ _ \ '__/ _ \ #
# | \__/\ (_) | | | | | | | | | | | (_| | | | | (_| \__ \ /\__/ / || (_| | |  | |_  | | | |  __/ | |  __/ #
#  \____/\___/|_| |_| |_|_| |_| |_|\__,_|_| |_|\__,_|___/ \____/ \__\__,_|_|   \__| \_| |_/\___|_|  \___| #
#                                                                                                         #
# ------------------------------------------------------------------------------------------------------- #

### Selecting which calculation types will be performed ###
calculations = {'undeformed': '', 'Applied_Field': '', 'uniaxial': ['0.9', '0.95', '1.05', '1.1']}

### Setting which database we work with ###
wdatadir = 'D:/MCES/TESTS/Gilles/datadir/'
wdatalist= 'D:/MCES/TESTS/Gilles/datalist_updated_no.duplicates_sieved.mag.sites_sieved.mag.active_sieved.mag.field.csv'
wdatalist = 'D:/MCES/TESTS/Gilles/datalist_updated_no.duplicates_sieved.mag.sites_sieved.mag.active_sieved.mag.field_beforeRunning.csv'
vasp_results_dir = 'D:/MCES/TESTS/Gilles/donedir'

# status_before(wdatalist_before, calculations)

# wdatalist = 'D:/MCES/TESTS/Gilles/datalist_updated_no.duplicates_sieved.mag.sites_sieved.mag.active_sieved.mag.field_beforeRunning.csv'
# status_after(wdatalist, vasp_results_dir, calculations)

wdatalist = 'D:/MCES/TESTS/Gilles/datalist_updated_no.duplicates_sieved.mag.sites_sieved.mag.active_sieved.mag.field_beforeRunning_afterRun.csv'
sieve_for_success(wdatalist, calculations)
