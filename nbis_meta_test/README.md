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
- `annotate` to annotate open reading frames called on assembled contigs using settings defined in the config file. It also quantifies genes and features, producing normalized and raw counts. When specified `True` in taxonomy, this step would include the taxonomy analysis.
- run statistical analysis on the outputs from `annotate` (including quantification and taxonomy reports).


### 1. QC 
To run the QC process, the command was used as: (with the config file `config.yaml`)
```shell
cd nbis-meta
conda activate nbis-meta
# the temp directory should be created mannually
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
Since we specified `True` in taxonomy analysis under the command `annotate`, the command would also run the taxonomy step automatically. To ensure the snakemake workflow could run successfully, we chose to download the UniRef100 database manually under the target directory.
```shell
wget https://ftp.uniprot.org/pub/databases/uniprot/uniref/uniref100/uniref100.fasta.gz /resources/uniref100/
```


### 4. Statistical Analysis



