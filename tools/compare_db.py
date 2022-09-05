import pandas as pd
import tqdm


def intersection(datalist1, datalist2):
    """Checks how many entries are present in both databases based on the field 'compound'"""
    df1 = pd.read_csv(datalist1, index_col=0, sep=',')
    df2 = pd.read_csv(datalist2, index_col=0, sep=',')
    s1 = pd.merge(df1, df2, how='inner', on=['compound'])

    s1.to_csv('X:/MCES/compare/diff.csv')
    print(df1.shape[0])
    print(df2.shape[0])
    print('intersection is ', s1.shape[0])


def apply_column(datalist1, datalist2, column):
    """Adds the values of chosen column from datalist2 to datalist1.
     Useful for comparing before/after of some parameter"""
    df1 = pd.read_csv(datalist1, index_col=0, sep=',')
    df2 = pd.read_csv(datalist2, index_col=0, sep=',')

    with tqdm.tqdm(total=len(df1.index)) as pbar:        # A wrapper that creates nice progress bar
        pbar.set_description("Processing datalist")
        for item in df1.index.tolist():
            try:
                pbar.update(1)  # Updating progress bar at each step
                df1.loc[item, column + '_' + datalist2] = df2.loc[item, column]
            except: pass

    df1.to_csv(datalist1.replace(".csv", '_' + column + '.compare' + '.csv'))

# ------------------------------------------------------------------------------------------------------- #
#  _____                                           _       _____ _             _     _   _                #
# /  __ \                                         | |     /  ___| |           | |   | | | |               #
# | /  \/ ___  _ __ ___  _ __ ___   __ _ _ __   __| |___  \ `--.| |_ __ _ _ __| |_  | |_| | ___ _ __ ___  #
# | |    / _ \| '_ ` _ \| '_ ` _ \ / _` | '_ \ / _` / __|  `--. \ __/ _` | '__| __| |  _  |/ _ \ '__/ _ \ #
# | \__/\ (_) | | | | | | | | | | | (_| | | | | (_| \__ \ /\__/ / || (_| | |  | |_  | | | |  __/ | |  __/ #
#  \____/\___/|_| |_| |_|_| |_| |_|\__,_|_| |_|\__,_|___/ \____/ \__\__,_|_|   \__| \_| |_/\___|_|  \___| #
#                                                                                                         #
# ------------------------------------------------------------------------------------------------------- #

dat1 = 'D:/MCES/COD/step2_beforeRunning_afterRun_success_sieved_out.csv'
dat2 = 'D:/MCES/COD/COD_before_error_fixes.csv'

# intersection(dat1, dat2)
apply_column(dat1, dat2, 'Mel_full')
