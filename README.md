# Mapping the Cellular Etiology of Psychiatric Disorders

## Software Requirements
1. Python libraries: `h5py`, `numexpr`
2. R libraries: `tidyverse`, `rhdf5`, `AnnotationDbi`, `org.Hs.eg.db`, `dplyr`, `readr` 
3. [MAGMA v1.10](https://ctg.cncr.nl/software/magma)

## Data Requirements
1. [Siletti et al.'s single-cell RNAseq dataset](https://github.com/linnarsson-lab/adult-human-brain)
2. [MAGMA auxiliary files](https://ctg.cncr.nl/software/magma)
3. GWAS summary statistics

## Get MAGMA Inputs
Follow the following steps:
1. [Get the ln(1+x)-transformed cluster-by-gene matrix.](Preprocessing_Siletti/create_matrices/Siletti_create_L2-log_dataset.py)
2. [Preprocess the matrix and calculate specificity.](Preprocessing_Siletti/create_magma_inputs/get_Siletti_continuous_input.md)

## Run MAGMA
First, confirm the following items:
1. The summary statistics are from a single population that matches MAGMA's auxiliary data.
2. The summary statistics are of the same genome build as MAGMA's auxiliary files.
3. If the summary statistics do not contain a SNP ID column, obtain the SNP IDs from the chromosomal and base pair positions using a reference file of the same genome build.

Then, follow the steps below to run MAGMA (scripts to be modified accordingly):
1. [Annotate and conduct a gene analysis.](MAGMA/1.annotationAndGeneAnalysis.sh)
2. [Run a gene property analysis.](MAGMA/2.genePropertyAnalysis.sh)

### Conditional Analysis
To run a pairwise conditional analysis on clusters after the steps above, follow these steps (scripts to be modified accordingly):
1. To limit the computation time, [create a MAGMA input file with only top clusters](MAGMA/3.create_top_results_matrix.md), as indicated by the results from the previous step.
2. [Run a pairwise conditional analysis.](MAGMA/4.conditionalAnalysis.sh)
3. [Conduct a forward stepwise selection](MAGMA/5.forward_selection_condition_results.md) to arrive at a list of independent clusters.
