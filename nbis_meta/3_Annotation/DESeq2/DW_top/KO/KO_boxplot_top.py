# -*- coding: utf-8 -*-
#!/usr/bin/python3
"""

Title: KO_boxplot_top.py
Date:2023-05-26
Author: Guyuan Tang

Description: The program was designed to select interested EC number and generate a box plot to show the differences of the distribution of the mean normalized counts between dry and wet samples.

List of functions:
    boxplot_KO

Procedure:
    1. According to the EC number entered by the user, select the normalized counts from the dataset.
    2. Create a new dataframe with the normalized counts and their sample conditions (dry / wet).
    3. Draw the boxplot based on the dataframe.
    

Usage: python KO_boxplot_top.py 1.14.13.25


"""

import matplotlib.pyplot as plt
import re
import pandas as pd
import numpy as np
import seaborn as sns
import sys

# define a function to generate boxplot
def boxplot_KO(KO_dict):
    for KO in KO_dict.keys():
        counts = KO_dict[KO]
        conditions = ["Dry","Dry","Dry","Wet","Wet","Wet","Wet","Wet"]
        a = {"counts":counts, "conditions":conditions}
        df = pd.DataFrame(a)
        out_name = KO + ".png"
        plt.figure()
        sns.boxplot(data=df, x="conditions", y="counts").set(title = KO)
        plt.savefig(out_name, dpi = 300)

# Main code
## laod the data
sig_KO_tab = "diffAbundanceSig_top.tab"
KO_counts_tab = "ko.normalised.counts.tab"
queryEC = sys.argv[1]

with open(sig_KO_tab, 'r') as sig_KO_tab, open(KO_counts_tab, 'r') as KO_counts_tab:
    # create an empty dictionary to store the information
    KO_dict = {}
    # define the regrex pattern
    queryPa = ''
    for i in range(len(queryEC)):
        if queryEC[i] == '.':
            queryPa += '\.'
        else:
            queryPa += queryEC[i]
    pattern = queryPa + '[^\d]'
    # search for the EC number in the sig_KO_tab to get the KOs
    for line in sig_KO_tab:
        line = line.strip().split(sep="\t")[0]
        if re.search(pattern, line):
            KO = line[:6]
            # search for normalized counts in KO_counts_tab
            for content in KO_counts_tab:
                content = content.strip()
                if content.startswith(KO):
                    count_info = content.split(sep="\t")[1:]
                    count_info = list(map(float, count_info))
                    KO_dict[KO] = count_info
                    
    # draw the boxplot
    boxplot_KO(KO_dict)
    
        
    
