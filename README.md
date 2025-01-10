# Mapping the Cellular Etiology of Schizophrenia and Diverse Brain Phenotypes
This is the code for the following publication:
Laramie E Duncan*, Tayden Li*, Madeleine Salem, Will Li, Leili Mortazavi, Hazal Senturk, Naghmeh Shargh, Sam Vesuna, Hanyang Shen, Jong Yoon, Gordon Wang, Jacob Ballon, Longzhi Tan, Brandon Scott Pruett, Brian Knutson, Karl Deisseroth, William J Giardino. Mapping the Cellular Etiology of Schizophrenia and Diverse Brain Phenotypes (in revision).


## Software Requirements
1. Python libraries: `h5py`, `numexpr`
2. R libraries: `tidyverse`, `rhdf5`, `AnnotationDbi`, `org.Hs.eg.db`, `dplyr`, `readr` 
3. [MAGMA v1.10](https://cncr.nl/research/magma/)

## Data Requirements
1. GWAS summary statistics
   - The schizophrenia summary statistics used in the example code can be downloaded from the [PGC website](https://figshare.com/articles/dataset/scz2022/19426775)
   - File name: `PGC3_SCZ_wave3.european.autosome.public.v3.vcf.tsv`
   - File md5sum = `2955c20b93f62607d650c83e7c41b0c7`
2. [Siletti et al.'s single-cell RNAseq dataset](https://github.com/linnarsson-lab/adult-human-brain)
3. [MAGMA auxiliary files](https://cncr.nl/research/magma/) of the same genome build and ancestry as GWAS summary statistics to be used
    - The schizophrenia summary statistics used in the example code are build 37 and only include individuals of European ancestry
    - [Gene locations, build 37 file](https://vu.data.surfsara.nl/index.php/s/Pj2orwuF2JYyKxq): `NCBI37.3.gene.loc`, md5sum = `7f6ffd6abb0b02c94a629e74c036d5f6`
    - [European ancestry reference data](https://vu.data.surfsara.nl/index.php/s/VZNByNwpD8qqINe): download the folder `g1000_eur`, which contains
      - `g1000_eur.bed`, md5sum = `543f92fba6ffe987dabf848b169c968d`
      - `g1000_eur.bim`, md5sum = `48fbc06e2ec736fbcc58989cb1727322`
      - `g1000_eur.fam`, md5sum = `6935b96c7e6f2a8f6122a049969bd0b8`
      - `g1000_eur.synonyms`, md5sum = `2807fd31bf0e92685c63629dc49b6574`
    - Put these 5 files in a directory called `aux` to use when running MAGMA

## Get MAGMA Inputs
Follow the following steps:
1. [Get the ln(1+x)-transformed cluster-by-gene matrix.](Preprocessing_Siletti/create_matrices/Siletti_create_L2-log_dataset.py)
2. [Preprocess the matrix and calculate specificity.](Preprocessing_Siletti/create_magma_inputs/get_Siletti_continuous_input.md)

## Run MAGMA
First, confirm the following items:
1. The summary statistics are from a single population that matches MAGMA's auxiliary data.
2. The summary statistics are the same genome build as MAGMA's auxiliary files.
3. If the summary statistics do not contain a SNP ID column, obtain the SNP IDs from the chromosomal and base pair positions using a reference file of the same genome build.

Then, follow the steps below to run MAGMA (scripts to be modified accordingly):
1. Create a SNP location file (`snploc_{GWAS_file_name}`) that contains three columns of the GWAS summary statistics in the following order: SNP ID, chromosome, and base pair position.
2. Annotate and conduct a gene analysis.
     Example code is provided [here](MAGMA/1.annotationAndGeneAnalysis.sh). The annotation step requires SNP location files created earlier, while the gene analysis step requires original GWAS files. Please refer to the MAGMA manual for different specification options for sample size and more.
4. Run a gene property analysis.
     Example code is provided [here](MAGMA/2.genePropertyAnalysis.sh). Please also refer to the MAGMA manual for the usage of different flag options.

### Conditional Analysis
To run a pairwise conditional analysis on clusters after the steps above, follow these steps (scripts to be modified accordingly):
1. To limit the computation time, [create a MAGMA input file with only top clusters](MAGMA/3.create_top_results_matrix.md), as indicated by the results from the previous step.
2. [Run a pairwise conditional analysis.](MAGMA/4.conditionalAnalysis.sh)
3. [Conduct a forward stepwise selection](MAGMA/5.forward_selection_condition_results.md) to arrive at a list of independent clusters.
