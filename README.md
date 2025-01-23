# Mapping the Cellular Etiology of Schizophrenia and Diverse Brain Phenotypes
Laramie E Duncan*, Tayden Li*, Madeleine Salem, Will Li, Leili Mortazavi, Hazal Senturk, Naghmeh Shargh, Sam Vesuna, Hanyang Shen, Jong Yoon, Gordon Wang, Jacob Ballon, Longzhi Tan, Brandon Scott Pruett, Brian Knutson, Karl Deisseroth, William J Giardino. 

## Summary 
This repository presents a comprehensive data-driven approach to unraveling the cellular and molecular underpinnings of psychiatric disorders, with a focus on schizophrenia. By integrating single nuclei RNA sequencing (snRNAseq) data with genome-wide association studies (GWAS) results, our work identifies and characterizes specific brain cell types implicated in disease etiology. Our findings not only validate previously reported associations but also reveal novel cell type contributions, achieving greater molecular specificity than earlier studies. This project further lays the foundation for a cell-type based classification system and offers a strategic roadmap for drug repurposing, novel therapeutic development, and personalized treatment strategies for psychiatric and other complex brain disorders.

<img width="822" alt="Screenshot 2025-01-13 at 10 59 25 AM" src="https://github.com/user-attachments/assets/e06c8a9d-cedd-475b-ab2e-8436ec598138" />


> [!WARNING]
> In our analysis, the term "cluster" as used in the Siletti 2023 paper is referred to as "cell types". Please keep this terminology difference in mind when interpreting the results and comparing findings with the original publication.

## Environment Setup
Before running the analysis, please ensure that you have the following software and libraries installed.
1. Python libraries: `h5py`, `numexpr`
2. R libraries: `tidyverse`, `rhdf5`, `AnnotationDbi`, `org.Hs.eg.db`, `dplyr`, `readr` 
3. [MAGMA v1.10](https://cncr.nl/research/magma/)

## Data Requirements
1. GWAS summary statistics
   - The schizophrenia summary statistics used in the example code can be downloaded from the [PGC website](https://figshare.com/articles/dataset/scz2022/19426775)
   - File name: `PGC3_SCZ_wave3.european.autosome.public.v3.vcf.tsv`
   - File md5sum = `2955c20b93f62607d650c83e7c41b0c7`
2. [Siletti et al.'s single-cell RNAseq dataset](https://github.com/linnarsson-lab/adult-human-brain)
   - File name: `adult_human_20221007.loom`
   * Note: this is a large file and can take ~1 hour to download.
3. [MAGMA auxiliary files](https://cncr.nl/research/magma/) of the same genome build and ancestry as GWAS summary statistics to be used
   - The schizophrenia summary statistics used in the example code are build 37 and only include individuals of European ancestry
   - [Gene locations, build 37 file](https://vu.data.surfsara.nl/index.php/s/Pj2orwuF2JYyKxq): `NCBI37.3.gene.loc`
   - [European ancestry reference data](https://vu.data.surfsara.nl/index.php/s/VZNByNwpD8qqINe): download the folder `g1000_eur`, which contains
      - `g1000_eur.bed`
      - `g1000_eur.bim`
      - `g1000_eur.fam`
      - `g1000_eur.synonyms`
   - Put these 5 files in a directory called `aux` to use when running MAGMA
4. Pre-processing auxiliary files
   - For step 2 of "Get MAGMA Inputs" (Preprocess the matrix and calculate specificity), remove the extended MHC (chromosome 6, 25Mb to 34Mb) from the auxiliary gene locations file
   - For build 37, you can download the file `NCBI37.3.gene.loc.extendedMHCexcluded` [here](https://github.com/jbryois/scRNA_disease/blob/master/Code_Paper/Data/NCBI/NCBI37.3.gene.loc.extendedMHCexcluded) from Bryois et al.
  
> [!IMPORTANT]  
> Ensure Data Compatibility. Please confirm the following items before proceeding to the MAGMA analysis:
> 1. The summary statistics are from a non-admixed population that matches MAGMA's auxiliary data.
> 2. The summary statistics are the same genome build as MAGMA's auxiliary files.

## Get MAGMA Inputs
Follow these steps:
1. [Get the ln(1+x)-transformed cluster-by-gene matrix.](Preprocessing_Siletti/create_matrices/Siletti_create_L2-log_dataset.py)
2. [Preprocess the matrix and calculate specificity.](Preprocessing_Siletti/create_magma_inputs/get_Siletti_continuous_input.md)
3. Check that your GWAS summary statistics file has the following columns: `rsID`, `chromosome`, `base pair location`, `p-value`, and `n`
   - `n` = the number of cases + controls for the phenotype and ancestry of interest
   - If your file does not contain a `rsID` column, obtain the rsID numbers from the chromosomal and base pair positions using a reference file of the same genome build
   - Include a header row on this file
4. Create a SNP location file (`snploc_{GWAS_file_name}`)
   - This file should contain three columns of the GWAS summary statistics in the following order: `rsID`, `chromosome`, and `base pair position`
   - Remove any header rows from this file

## Run MAGMA
1. Annotate and conduct a gene analysis.
     Example code is provided [here](MAGMA/1.annotationAndGeneAnalysis.sh). The annotation step requires SNP location files created earlier, while the gene analysis step requires original GWAS files. Please refer to the MAGMA manual for different specification options for sample size and more.
2. Run a gene property analysis.
     Example code is provided [here](MAGMA/2.genePropertyAnalysis.sh). Please also refer to the MAGMA manual for the usage of different flag options.

## Conditional Analysis
To run a pairwise conditional analysis on clusters after the steps above, follow these steps (scripts to be modified accordingly):
1. To limit the computation time, [create a MAGMA input file with only top clusters](MAGMA/3.create_top_results_matrix.md), as indicated by the results from the previous step.
2. [Run a pairwise conditional analysis.](MAGMA/4.conditionalAnalysis.sh)
3. [Conduct a forward stepwise selection](MAGMA/5.forward_selection_condition_results.md) to arrive at a list of independent clusters.

## Cite us
```
@article {Duncan2025,
	title = {Mapping the Cellular Etiology of Schizophrenia and Diverse Brain Phenotypes},
	author = {Duncan, Laramie E and Li, Tayden and Salem, Madeleine and Li, Will and Mortazavi, Leili and Senturk, Hazal and Shargh, Naghmeh and Vesuna, Sam and Shen, Hanyang and Yoon, Jong and Wang, Gordon and Ballon, Jacob and Tan, Longzhi and Pruett, Brandon Scott and Knutson, Brian and Deisseroth, Karl and Giardino, William J},
	year = {2025},
	journal = {Nature Neuroscience},
	doi={10.1038/s41593-024-01834-w},
	url={https://www.nature.com/articles/s41593-024-01834-w}
}
```
