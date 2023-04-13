# README for Testing NBIS-meta Workflow
Author: Guyuan Tang  
Date: 13 Apr 2023

## Description
Before running nbis-meta on the required samples under the project, we tested the workflow to manage the settings and the environment.  
Samples and the config.yaml for this Snakemake were provided in the folder `SampleList`. The `config.yaml` was used throughout the testing process.

## Steps
- `qc` to run the pre-processing and to generate the quality report on the samples.
- `assemble` to run the metagenomic assembly on the preprocessed samples.

### 1. QC
 To run the QC process, the command was used as:
```shell
cd Project37
cd nbis-meta
# the temp directory should be created mannually
mkdir test_temp
snakemake --use-conda --configfile ./config/config.yaml --cores 60 qc
```
The final multiqc report `samples_report.html` was provided under the folder `1_QC`, including 27 samples.

### 2. Assembly
To run the assembly process, the command was used as:
```shell
snakemake --use-conda --configfile ./config/config.yaml --cores 60 assemble
```


