# Title: stat_KO_DW.R
# Author: Guyuan Tang
# Date: 2023-05-11

# Description: using DESeq2 to analyse the raw counts of KO annotations of different metagenome samples.
# Packages: DESeq2, gplots, RColorBrewer, dplyr, vegan

## Load the packages
library(DESeq2)
library(gplots)
library(RColorBrewer)
library(dplyr)
library(vegan)
library(ggplot2)

## Load the data
count_mat = as.matrix(read.csv('../kos.parsed.counts.tsv', sep = '\t', row.names = 1))
# remove the samples that are not required here
# here we remove the mesic IDs and IDs without metadata
# and keep the IDs collected at bottom
count_mat = subset(count_mat, select = c(ID18_1, ID12_1, ID20_1, ID34_1, ID22_1, ID36_1, ID32_1, ID38_1))
mode(count_mat) <- "integer"   # Convert to integer


## Create the metadata table
sampleNames = c("ID18_1", "ID12_1", "ID20_1", "ID34_1", "ID22_1", "ID36_1", "ID32_1", "ID38_1")
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
