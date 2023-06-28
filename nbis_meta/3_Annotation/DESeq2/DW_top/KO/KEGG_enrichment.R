# Title: KEGG_enrichment.R
# Author: Guyuan Tang
# Date: 2023-05-18

# Description: the program was designed to conduct the pathway enrichment analysis on the significantly different KO list.

## Install the package clusterProfiler
BiocManager::install("clusterProfiler")
library(clusterProfiler)
library(enrichplot)

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
enrichplot::emapplot(KEGG2,showCategory =50, color = "p.adjust", layout = "kk")
dev.off()
