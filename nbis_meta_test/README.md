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
- `assemble` to run the metagenomic assembly on the preprocessed samples.

### 1. QC
Since the pre-processing step would be the same for both assemblers, we used the config file `config.yaml` here.  
To run the QC process, the command was used as:
```shell
cd nbis-meta
conda activate nbis-meta
# the temp directory should be created mannually
mkdir test_temp
snakemake --use-conda --configfile ./config/config.yaml --cores 60 qc
```
The final multiqc report `samples_report.html` was provided under the folder `1_QC`, including 27 samples.

### 2. Assembly
The metagenomic assembly for both assemblers would be on each ID only as well as all the SF samples, which were specified in the `Abisko_SF2.tsv`.
#### (1) Using megahit as the assembler
To run the assembly process in megahit, the command was used as:
```shell
cd nbis-meta
conda activate nbis-meta
snakemake --use-conda --configfile ./config/config.yaml --cores 60 assemble
```
##### using metaquast to check the quality of assembly made by megahit
``` shell
cd test_results/assembly/
metaquast ./ID11/final_contigs.fa ./ID12/final_contigs.fa ./ID13/final_contigs.fa ./ID14/final_contigs.fa ./ID15/final_contigs.fa ./ID17/final_contigs.fa ./ID18/final_contigs.fa ./ID19/final_contigs.fa ./ID20/final_contigs.fa ./ID21/final_contigs.fa ./ID22/final_contigs.fa ./ID23/final_contigs.fa ./ID24/final_contigs.fa ./ID25/final_contigs.fa ./ID26/final_contigs.fa ./ID27/final_contigs.fa ./ID28/final_contigs.fa ./ID29/final_contigs.fa ./ID30/final_contigs.fa ./ID31/final_contigs.fa ./ID32/final_contigs.fa ./ID33/final_contigs.fa ./ID34/final_contigs.fa ./ID35/final_contigs.fa ./ID36/final_contigs.fa ./ID37/final_contigs.fa ./ID38/final_contigs.fa ./SF/final_contigs.fa -t 40 -o ../quast_results/
```

#### (2) Using metaspades as the assembler
To run the assembly process in metaspades, the command was used as:
```shell
cd nbis-meta
conda activate nbis-meta
# create the temporary directory manually
mkdir spades_test_temp
snakemake --use-conda --configfile ./config/config_metaspades.yaml --cores 60 assemble
```

