# README for Testing NBIS-meta Workflow
Author: Guyuan Tang  
Date: 13 Apr 2023

## Description
Before running nbis-meta on the required samples under the project, we tested the workflow to manage the settings and the environment.  
Samples and the config.yaml for this Snakemake were provided in the folder `SampleList`. The `config.yaml` and `config_metaspades.yaml` were used throughout the testing process with only difference in the choice of assemblers.

## Installation of NBIS-meta
The detailed installation guideline could be found in the github of NBIS-meta: https://github.com/NBISweden/nbis-meta.
```shell
cd Project37
git clone https://github.com/NBISweden/nbis-meta.git
conda env create -f environment.yml
conda activate nbis-meta
```


## Steps
- `qc` to run the pre-processing and to generate the quality report on the samples.
- `assemble` to run the metagenomic assembly (megahit) on the preprocessed samples. `metaquast` to compare the assembly results.
- `annotate` to annotate open reading frames called on assembled contigs using settings defined in the config file. It also quantifies genes and features, producing normalized and raw counts. (When specified `True` in taxonomy, this step would include the taxonomy analysis.)
- `taxonomy` to assign taxonomy to assembled contigs and ORFs called on these contigs.
- run statistical analysis on the outputs from `annotate` (including quantification reports) to investigate the differences on pathway enrichment between dry and wet samples from top and bottom layers separately, in order to discover the differences in functional compositions. The program languages included Python and R.
- run statistical analysis on the outputs form `taxonomy` to discover the differences in community compositions between dry and wet samples from top and bottom layers separately.


### 1. QC 
To run the QC process, the command was used as: (with the config file `config.yaml`)
```shell
cd nbis-meta
conda activate nbis-meta
# the temp directory should be created manually
mkdir test_temp
snakemake --use-conda --configfile ./config/config.yaml --cores 60 qc
```
The final multiqc report `samples_report.html` was provided under the folder `1_QC`, including 27 samples.

### 2. Assembly
The metagenomic assembly would be on each ID only as well as all the SF samples, which were specified in the `Abisko_SF1.tsv`. 

#### Using megahit as the assembler
To run the assembly process in megahit, the command was used as:
```shell
cd nbis-meta
conda activate nbis-meta
snakemake --use-conda --configfile ./config/config.yaml --cores 60 assemble
```
Megahit used the `Abisko_SF1.tsv` as the assemble sample setting.

##### Using metaquast to check the quality of assembly made by megahit
The comparison was made between the assembly on each ID and the assembly on all of the IDs, without reference. The final report `report.html` was saved to `2_Assembly_megahit/metaquast`.
``` shell
conda activate quast
cd test_results/assembly/
metaquast ./ID11/final_contigs.fa ./ID12/final_contigs.fa ./ID13/final_contigs.fa ./ID14/final_contigs.fa ./ID15/final_contigs.fa ./ID17/final_contigs.fa ./ID18/final_contigs.fa ./ID19/final_contigs.fa ./ID20/final_contigs.fa ./ID21/final_contigs.fa ./ID22/final_contigs.fa ./ID23/final_contigs.fa ./ID24/final_contigs.fa ./ID25/final_contigs.fa ./ID26/final_contigs.fa ./ID27/final_contigs.fa ./ID28/final_contigs.fa ./ID29/final_contigs.fa ./ID30/final_contigs.fa ./ID31/final_contigs.fa ./ID32/final_contigs.fa ./ID33/final_contigs.fa ./ID34/final_contigs.fa ./ID35/final_contigs.fa ./ID36/final_contigs.fa ./ID37/final_contigs.fa ./ID38/final_contigs.fa ./SF/final_contigs.fa -t 40 -o ../quast_results/
```


### 3. Annotation
While running the annotation on all of the samples but not on each ID, the samples were specified in the `Abisko_SF2.tsv` instead.  
The annotation step included the `quantification` step automatically.  
To run the annotation step, the following commands were used:
```shell
conda activate nbis-meta
cd nbis-meta
snakemake --use-conda --configfile ./config/config.yaml --cores 60 annotate
```
However, during the process, when producing the output file `modules.parsed.counts.tsv`, there were duplicated row names ("Phthalate degration, phthatlate => protocatechuate [Path:map00624 map01220 map01100 map01120]") which caused errors because the R script `edger.R` within the snakemake did not accept the duplicated row names.  
  
And thus, we mannually added "1" and "2" at the ends of these two row names to avoid the error messages. As the result, in `modules.parsed.counts.tsv`, row name for 255 would be "Phthalate degration, phthatlate => protocatechuate [Path:map00624 map01220 map01100 map01120]1". And row name for 256 would be "Phthalate degration, phthatlate => protocatechuate [Path:map00624 map01220 map01100 map01120]2".  
  
Then we continued the nbis-meta snakemake again after manually correting the duplicated row names.
```shell
snakemake --unlock
snakemake --rerun-incomplete --use-conda --configfile ./config/config.yaml --cores 60 annotate
```
If we specified `True` in taxonomy analysis under the command `annotate`, the command would also run the taxonomy step automatically. 


### 4. Taxonomy
However, to ensure the snakemake workflow could run successfully, we chose to download the UniRef100 database manually under the target directory.
```shell
cd nbis-meta
wget https://ftp.uniprot.org/pub/databases/uniprot/uniref/uniref100/uniref100.fasta.gz /resources/uniref100/
```
Before running the `taxonomy`, in order to improve the efficiency and speed of the programme, we altered the number of "threads" in the snakemake file `taxonomy.smk`, which were set to **200**. Then we ran the `taxonomy` command as the following:
```shell
snakemake --use-conda configfile ./config/config.yaml --cores 200 taxonomy
```


### 5. Statistical Analysis on Annotation Outputs
To analyze the differences in functional compositions between dry and wet samples from top and bottom layers separately (here the information on KO was chosen), we wrote programs to:  
- `annotation_stats.py` to generate basic statistics of the annotation outputs, including the mean counts per row, number of rows, percentage of sparse in the data, and the sum of per row as well as per ID;
- `stat_KO_DW.R` used DESeq2 to analyse the raw counts of KO annotations of different metagenome samples, including calculating the normalized counts and creating the PCA and the NMDS plots;
- `KO_Parser.py` to select and calculate the mean of significant KO counts of selected group. With the list of significantly abundant KO and their means, the program would assign colours to each level of counts for later KO mapper (https://www.genome.jp/kegg/mapper/color.html);
- `KEGG_enrichment.R` to conduct the pathway enrichment analysis on the significantly different KO list;
- `KO_boxplot_top.py` and `KO_boxplot_bottom.py` to select interested EC number and generate a box plot to show the differences of the distribution of the mean normalized counts between dry and wet samples.  

#### Information on `annotation_stats.py`
We used the following command to run the program:
```shell
python annotation_stats.py kos.parsed.TMM.tsv
```
and get the output file `kos.parsed.TMM.txt` with the basic information on the annotation outputs.  
The script was as the followings:
```python
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
```

#### Information on `stat_KO_DW.R`
The script was used in the directory that would be used to store the results (including the normalized counts table, PCA plot, and NMDS plot) from DESeq2 (v1.40.1).
```R
## Load the packages
library(DESeq2)
library(gplots)
library(RColorBrewer)
library(dplyr)
library(vegan)
library(ggplot2)

## Load the data
count_mat = as.matrix(read.csv('../kos.parsed.counts.tsv', sep = '\t', row.names = 1)) # the directory could be altered 
# remove the samples that are not required here
# here we remove the mesic IDs and IDs without metadata
# and keep the IDs collected at top
count_mat = subset(count_mat, select = c(ID11_1, ID17_1, ID19_1, ID35_1, ID21_1, ID33_1, ID37_1, ID31_1))
# if the samples were from the bottom
#count_mat = subset(count_mat, select = c(ID18_1, ID12_1, ID20_1, ID34_1, ID22_1, ID36_1, ID32_1, ID38_1))
mode(count_mat) <- "integer"   # Convert to integer


## Create the metadata table
sampleNames = c("ID11_1", "ID17_1", "ID19_1", "ID35_1", "ID21_1", "ID33_1", "ID37_1", "ID31_1")
# if the samples were from the bottom
#sampleNames = c("ID18_1", "ID12_1", "ID20_1", "ID34_1", "ID22_1", "ID36_1", "ID32_1", "ID38_1")
sampleConditions = c("Dry", "Dry", "Dry", "Wet", "Wet", "Wet", "Wet", "Wet")
sampleTable <- data.frame(condition = as.factor(sampleConditions))
row.names(sampleTable) <- sampleNames # Set the first column as row names


## Create the DESeq object
dds = DESeqDataSetFromMatrix(countData = count_mat, colData = sampleTable, design = ~ condition)

## Normalized the counts within the DESeq object
dds = estimateSizeFactors(dds)
sizeFactors(dds)
normalized_counts = as.data.frame(counts(dds, normalized = TRUE))
summary(normalized_counts)
# save the normalized table
nor_names = rownames(normalized_counts)
normalized_out = cbind(nor_names, data.frame(normalized_counts, row.names = NULL))
write.table(normalized_out, file = 'ko.normalised.counts.tab', sep = '\t', quote = FALSE, row.names = FALSE)

## PCA
# regularized log transform
rld = rlogTransformation(dds, blind = TRUE)
# PCA plot
pdf('fig1.PCA.pdf')
print(plotPCA(rld, intgroup = c('condition')))
dev.off()

## NMDS
normalized_counts = as.matrix(normalized_counts)
nor_count_t = t(normalized_counts)
# calculate relative abundance
count_rel = decostand(nor_count_t, method = 'total')
# calculate distance matrix
dist_mat = vegdist(count_rel, method = 'bray')
dist_mat = as.matrix(dist_mat, labels = TRUE)
# run the NMDS
nmds = metaMDS(dist_mat, distance = 'bray')
# extract NMDS scores
data.scores = as.data.frame(scores(nmds))
data.scores$condition = sampleTable$condition
# plot the NMDS results with ggplot2
xx = ggplot(data.scores, aes(x = NMDS1, y = NMDS2)) + 
  geom_point(size = 2, aes( shape = condition))+ 
  theme(axis.text.y = element_text(colour = "black", size = 10), 
        axis.text.x = element_text(colour = "black", size = 10), 
        legend.text = element_text(size = 10, colour ="black"), 
        legend.position = "right", axis.title.y = element_text(size = 11), 
        axis.title.x = element_text(size = 11, colour = "black"), 
        legend.title = element_text(size = 11, colour = "black"), 
        panel.background = element_blank(), panel.border = element_rect(colour = "black", fill = NA, size = 1.2),
        legend.key=element_blank()) + 
  labs(x = "NMDS1", y = "NMDS2", shape = "Condition")
ggsave("fig2.NMDS.png")



## Differential abundance
# run the analysis
dds = DESeq(dds)
# specify the contrast
contrast_pr <- c("condition", "Dry", "Wet")
# extract the results for the specified contrast
res_table <- results(dds, contrast=contrast_pr)
# sort results based on the padj
res_table <- res_table[order(res_table$padj),]
head(res_table)
# save the information of res_table
anno_names = rownames(res_table)
res_table = cbind(anno_names, data.frame(res_table, row.names = NULL))
write.table(res_table, file = 'diffAbundance.tab', sep = '\t', quote = FALSE, row.names = FALSE)
# only those with significant differences (padj<0.05)
resSig = subset(res_table, padj<0.05)
write.table(resSig, file = 'diffAbundanceSig.tab', sep = '\t', quote = FALSE, row.names = FALSE)
```

#### Information on `KO_Parser.py`
The program was used in the directory that stored the normalized counts table.  
For example, if we would like to see which KOs were significantly more abundant in wet samples from top (WT), we ran the following command:  
`python KO_Parser.py WT diffAbundanceSig_top.tab ko.normalised.counts.tab`  
(WT: wet from top, DT: dry from top, WB: wet from bottom, DB: dry from bottom)  
The outputs for this program would include two text files: the KOs and their mean in the selected group (e.g. `WT_KO_m.txt`) and the color list (e.g. `WT_color.txt`).  
  
The script was as the followings:
```python
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
            color = '#0076FF'
        elif 100 <= val < 500:
            color = '#9AC9FF'
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
```


#### Information on `KEGG_enrichment.R`
The scipt was ran in the directory that stored the outputs from `KO_Parser.py` (the KOs and their mean normalized counts, e.g. WT_KO_m.txt). The package clusterProfiler (v4.8.1) was used in this part.  
The outputs from this script would be an enrichment dotplot, an enrichment barplot, an enrich plot and an emapplot. We decided to use the dotplot and the barplot showing the significant enriched pathways.
```R
## Install the package clusterProfiler
BiocManager::install("clusterProfiler")
library(clusterProfiler)
library(enrichplot)

# the group (WT/DT/WB/DB) could be selected by the user.

## Retrieve the KO ids from the input file
KO_info <- read.csv('WT_KO_m.txt', header = FALSE, row.names = NULL, sep = "\t", col.names = c("KO","mean"))

## Enrichment analysis
res <- enrichKEGG(KO_info$KO, organism = "ko", pvalueCutoff = 0.05, pAdjustMethod = "BH", qvalueCutoff = 0.05)
# print the dot plot
pdf("WT_dotplot_bh.pdf")
dotplot(res)
dev.off()
# print the bar plot
pdf("WT_barplot_bh.pdf")
barplot(res,showCategory = 20,title = 'KEGG Pathway')
dev.off()
# print the enrich plot
pdf("WT_enrich.pdf")
enrichplot::cnetplot(res,circular=FALSE,colorEdge = TRUE)
dev.off()


KEGG2 <- pairwise_termsim(res)
# print the emapplot
pdf("WT_emap.pdf")
enrichplot::emapplot(KEGG2,showCategory = 50, color = "p.adjust", layout = "kk")
dev.off()
```

#### Information on `KO_boxplot_top.py` and `KO_boxplot_bottom.py`
These two scripts were mostly the same instead of changing the sample IDs. Therefore, the script here we used the top as the example, the output of which would be the boxplot on the selected EC number.  
The program should be used as: `python KO_boxplot_top.py 1.14.13.25`
```python
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
```

