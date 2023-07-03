# -*- coding: utf-8 -*-
#!/usr/bin/python3
"""

Title: KO_Parser.py
Date:2023-05-17
Author: Guyuan Tang

Description: The program was designed to select and calculate the mean of significant KO counts of selected group. With the list of significantly abundant KO and their means, the program would assign colours to each level of counts for later KO mapper.

List of functions:
    sample_mean, color_KO

Procedure:
    1. According to the group selected by the user, filter the significantly more abundant KO names into a list.
    2. Seach for the KO names in the normalised counts table to retrieve the sample counts, and calculate the mean. Print the KO and the mean counts to an output file for later application.
    3. Assign colours to each KO based on their counts.
    

Usage: python KO_Parser.py WT diffAbundanceSig_top.tab ko.normalised.counts.tab


"""

import sys
import pandas as pd

# Define a function to select the dry or wet samples from the counts table and calculate their means
def sample_mean(counts_tab, dw_group, KO_list):
    # create an empty dictionary to contain the KOs and their means
    KO_mean = {}
    # loop over the KO_list
    for KO_name in KO_list:
        if dw_group == 'D': # for dry samples
            df = counts_tab.iloc[:,0:3]
            df = df[df.index.str.startswith(KO_name)]
            m_val = round(float(df.mean(axis=1)),2)
            KO_mean[KO_name] = m_val
        else: # for wet samples
            df = counts_tab.iloc[:,3:]
            df = df[df.index.str.startswith(KO_name)]
            m_val = round(float(df.mean(axis=1)),2)
            KO_mean[KO_name] = m_val
    return KO_mean


# Define a function to assign colors to each KO
def color_KO(KO_mean):
    KO_color = {}
    # assign colors to different count category
    for KO, val in KO_mean.items():
        if 0 <= val < 100:
            color = '#9AC9FF'
        elif 100 <= val < 500:
            color = '#0076FF'
        elif 500 <= val < 1000:
            color = '#FAC2E3'
        elif 1000 <= val < 2000:
            color = '#FB7B9A'
        elif 2000 <= val < 3000:
            color = '#FF0000'
        elif val >= 3000:
            color = '#D31D1D'
        
        KO_color[KO] = color
        
    return KO_color


# Load the data
dw_group = sys.argv[1][0]
sig_KO = sys.argv[2]
counts_tab = pd.read_csv(sys.argv[3], sep='\t', header=0, index_col=0)
# Set the output file names
out_mean = sys.argv[1] + '_KO_m.txt'
out_color = sys.argv[1] + '_color.txt'

with open(sig_KO,'r') as sig_KO, open(out_mean, 'w') as out_mean, open(out_color, 'w') as out_color:
    # create an empty list to contain the name of filtered KO
    KO_list = []
    header = sig_KO.readline() # ignore the header line
    for line in sig_KO:
        line = line.strip().split("\t")
        KO_name = line[0][:6] # only remain the KO number
        lg2f = float(line[2])
        if dw_group == 'D': # for dry samples
            if lg2f >= 0:
                KO_list.append(KO_name)
        else: # for wet samples
            if lg2f < 0 :
                KO_list.append(KO_name)
    
    # calculate the mean
    KO_mean = sample_mean(counts_tab, dw_group, KO_list)
    # print the results
    for KO, m in KO_mean.items():
        print('{}\t{}'.format(KO,m), file=out_mean)
    
    # assign colors
    KO_color = color_KO(KO_mean)
    # print the results
    for KO, color in KO_color.items():
        print('{}\t{}'.format(KO, color), file=out_color)



