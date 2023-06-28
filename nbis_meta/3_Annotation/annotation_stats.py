#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Title: annotation_stats.py
Author: Guyuan Tang
Date: 2023-05-10

Description:
    The program was designed to generate basic statistics of the annotation outputs, including the mean counts per row, number of rows, percentage of sparse in the data, and the sum of per row as well as per ID. It could help the user to have a general sense of the count data and choose type of the dataset or the filters in the later research steps.

Procedure:
    1. read the tsv file into dataframe;
    2. generate the number of rows (annotations) and IDs;
    3. calculate the percentage of sparses in the data;
    4. calculate the mean counts per row;
    5. calculate the sum counts per row and per ID;
    6. print the results to the txt output file.

List of functions:
    None.

Usage:
    python annotation_stats.py input.tsv
"""

import pandas as pd
import numpy as np
import sys

# read the tsv annotation file into dataframe
in_file = sys.argv[1]
out_file = in_file[:-3] + 'txt'

anno_df = pd.read_csv(in_file, sep='\t', header=0, index_col=0)

with open(out_file, 'w') as out_file:
    # generate the number of rows and IDs
    num_row = len(anno_df)
    stat_df = anno_df.loc[:,anno_df.columns.str.startswith("ID")] # extract the sub-dataframe with only the counts on IDs
    num_ID = stat_df.shape[1]
    # print the results to output
    print('# number of annotations: {}\n# number of IDs: {}'.format(num_row, num_ID), file=out_file)

    # calculate the percentage of sparses in the data
    total_num = num_row * num_ID
    num_zero = 0
    for col_name in stat_df.columns:
        count = (stat_df[col_name] == 0).sum()
        num_zero += count
    percent_zero = round(num_zero / total_num * 100, 2)
    # print the result to output
    print('# %sparse: {}'.format(percent_zero), file=out_file)

    # set the dataframe printing setting
    pd.set_option('display.max_rows', None) # to print all the rows
    pd.set_option('display.float_format', lambda x: '%.2f' % x) # to print the readable numbers
    
    # calculate the mean and sum counts per row
    print('\n# the sum and mean counts per row', file=out_file)
    row_mean = stat_df.mean(axis=1)
    row_sum = stat_df.sum(axis=1)
    row_df = pd.concat([row_sum, row_mean], axis=1, ignore_index=False)
    row_df.columns = ['sum','mean']
    row_df = row_df.sort_values(by=['sum'], ascending=False)
    # print the results to output
    print(row_df, file = out_file)
    
    # calculate the sum and mean counts per ID
    print('\n# the sum counts per ID', file=out_file)
    col_mean = stat_df.mean(axis = 0)
    col_sum = stat_df.sum(axis = 0)
    col_df = pd.concat([col_sum, col_mean], axis=1, ignore_index=False)
    col_df.columns = ['sum','mean']
    col_df = col_df.sort_values(by=['sum'], ascending=False)
    # print the results to output
    print(col_df, file = out_file)
