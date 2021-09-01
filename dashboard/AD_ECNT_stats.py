#!/usr/bin/env python3

import argparse
import pandas as pd
import sys

def main():
    parser = argparse.ArgumentParser(description='Calculate distribution stats from step 9 output for AD and ECNT')
    parser.add_argument('step9_output', type=str)
    parser.add_argument('out_csv', type=str)
    args = parser.parse_args()

    print(args)

    step9_df = read_csv(args.step9_output)

    AD_metrics = get_AD_statistics(step9_df)
    ECNT_metrics = get_ECNT_statistics(step9_df)
    DP_metrics = get_DP_statistics(step9_df)
    out_df = pd.DataFrame({"AD": AD_metrics, 
        "ECNT" : ECNT_metrics,
        "DP" : DP_metrics})

    write_csv(out_df, args.out_csv)


# run pandas.read_csv with exception handling
def read_csv(path):
    try:
        f = pd.read_csv(path)
    except OSError:
        print("Could not open/read file:", fname)
        sys.exit(1)
    return f

# run pandas.DataFrame.to_csv with exception handling
def write_csv(df, path):
    try:
        f = open(path, 'w')
    except OSError:
        print("Could not open/read file:", path)
        sys.exit(1)
    df.to_csv(f)

#get summary statistics of ECNT values
def get_ECNT_statistics(df):
    return df['ECNT'].describe()

def get_from_stat1_stat2(df, metric):
    # have to look at the stat1 and stat2 columns
    # the stat2 column has different metrics separated by colons
    # the stat1 column tells us what order the metrics are in
    stat1_split = df['stat1'].str.split(":", expand=True)
    stat2_split = df['stat2'].str.split(":", expand=True)

    # we look at stat1 column to see where in stat2 we'll find the wanted metric
    wanted_vals = stat2_split[stat1_split == metric]
    # make sure there's only one AD for each row
    assert((wanted_vals.notnull().sum(axis='columns')==1).all())

    wanted_vals = wanted_vals.stack().dropna(how='all')
    return wanted_vals


# get summary statistics of AD values
def get_AD_statistics(df):
    # the AD metric are in the format "<num_reads_ref>,<num_reads_alt>" in stat2
    # get summary statistics of <num_reads_alt>.
    AD_vals = get_from_stat1_stat2(df, "AD")
    num_reads_alt = AD_vals.str.split(",", expand=True)[1].astype(int)
    return num_reads_alt.describe()

# get summary statistics of DP values
def get_DP_statistics(df):
    # the AD metric are in the format "<num_reads_ref>,<num_reads_alt>" in stat2
    # get summary statistics of <num_reads_alt>.
    DP_vals = get_from_stat1_stat2(df, "DP").astype(int)
    return DP_vals.describe()



if __name__ == "__main__":
    main()

