---
title: "Siletti 2022 scRNAseq Dataset - Get Specificity Matrix"
author: "Tayden Li"
output:
  html_document:
    keep_md: yes
---

Code is adapted from [Bryois (2020)](https://github.com/jbryois/scRNA_disease/tree/master).
The script produces continuous specificity MAGMA inputs, with specificity calculated on a subset of clusters.

# Load Data

The Siletti data was downloaded from [here](https://github.com/linnarsson-lab/adult-human-brain)

### Load necessary libraries


```r
library(tidyverse)
```

```
## ── Attaching packages ─────────────────────────────────────── tidyverse 1.3.2 ──
## ✔ ggplot2 3.4.0      ✔ purrr   0.3.5 
## ✔ tibble  3.1.8      ✔ dplyr   1.0.10
## ✔ tidyr   1.2.1      ✔ stringr 1.5.0 
## ✔ readr   2.1.3      ✔ forcats 0.5.2 
## ── Conflicts ────────────────────────────────────────── tidyverse_conflicts() ──
## ✖ dplyr::filter() masks stats::filter()
## ✖ dplyr::lag()    masks stats::lag()
```

```r
library("rhdf5")
```

### Load single cell dataset


```r
file="/oak/stanford/groups/laramied/LDSC/Bryois2020/Code_Paper/Data/Siletti/Siletti_L2-cluster-log1p_matrix.h5"
h5f <- H5Fopen(file)
exp <- as.data.frame(t(h5f$matrix))
colnames(exp) <- h5f$Cluster # change column name
exp$Gene <- sub("\\..*", "", h5f$Accession)
h5closeAll()
```

### Keep certain clusters only.


```r
exp <- exp %>% select(Cluster83:Cluster205, Gene)
```


Only keep genes with a unique name and tidy data.


```r
exp <- exp %>% add_count(Gene) %>% 
  filter(n==1) %>%
  select(-n) %>%
  gather(key = column,value=Expr,-Gene) %>%
  as.tibble()
```

```
## Warning: `as.tibble()` was deprecated in tibble 2.0.0.
## ℹ Please use `as_tibble()` instead.
## ℹ The signature and semantics have changed, see `?as_tibble`.
```

### Load gene coordinates

Load gene coordinates and extend upstream and downstream coordinates by 100kb.

File downloaded from MAGMA website (https://ctg.cncr.nl/software/magma).

Filtered to remove extended MHC (chr6, 25Mb to 34Mb).


```r
gene_coordinates <- 
  read_tsv("../Data/NCBI/NCBI37.3.gene.loc.extendedMHCexcluded",
           col_names = FALSE,col_types = 'cciicc') %>%
  mutate(start=ifelse(X3-100000<0,0,X3-100000),end=X4+100000) %>%
  select(X2,start,end,1) %>% 
  rename(chr="X2", ENTREZ="X1") %>% 
  mutate(chr=paste0("chr",chr))
```


### Get table for ENTREZ and ENSEMBL gene names.


```r
entrez_ensembl <- AnnotationDbi::toTable(org.Hs.eg.db::org.Hs.egENSEMBL)
```

```
## 
```

Only keep genes with a unique entrez and ensembl id.


```r
entrez_ensembl_unique_genes_entrez <- entrez_ensembl %>% count(gene_id) %>% filter(n==1)
entrez_ensembl_unique_genes_ens <- entrez_ensembl %>% count(ensembl_id) %>% filter(n==1)
entrez_ensembl <- filter(entrez_ensembl,gene_id%in%entrez_ensembl_unique_genes_entrez$gene_id & ensembl_id %in% entrez_ensembl_unique_genes_ens$ensembl_id)
colnames(entrez_ensembl)[1] <- "ENTREZ"
colnames(entrez_ensembl)[2] <- "Gene"
gene_coordinates <- inner_join(entrez_ensembl,gene_coordinates) %>% as.tibble()
```

```
## Joining, by = "ENTREZ"
```



# Transform Data

### Add cell cluster ID


```r
exp_lvl5 <- transform(exp,ClusterID=exp$column) %>%
  rename(Expr_sum_mean=Expr)
```

### Remove not expressed genes


```r
not_expressed <- exp_lvl5 %>% 
  group_by(Gene) %>% 
  summarise(total_sum=sum(Expr_sum_mean)) %>% 
  filter(total_sum==0) %>% 
  select(Gene) %>% unique() 

exp_lvl5 <- filter(exp_lvl5,!Gene%in%not_expressed$Gene)
```

### Scale data

Each cell type is scaled to the same total number of molecules. 


```r
exp_lvl5 <- exp_lvl5 %>% 
  group_by(ClusterID) %>% 
  mutate(Expr_sum_mean=Expr_sum_mean*1000/sum(Expr_sum_mean))
```


# Specificity Calculation

The specifitiy is defined as the proportion of total expression performed by the cell type of interest (x/sum(x)).

### cluster level


```r
exp_lvl5 <- exp_lvl5 %>% 
  group_by(Gene) %>% 
  mutate(specificity=Expr_sum_mean/sum(Expr_sum_mean)) %>% 
  ungroup()
```

### Get MAGMA genes

Only keep genes that are tested in MAGMA


```r
exp_lvl5 <- inner_join(exp_lvl5,gene_coordinates,by="Gene")
```

### Get number of genes

Get number of genes that represent 10% of the dataset


```r
n_genes <- length(unique(exp_lvl5$ENTREZ))
```

### Write MAGMA input file


```r
exp_conti_spe <- exp_lvl5 %>% select(ENTREZ, ClusterID, specificity) %>% spread(ClusterID, specificity)
colnames(exp_conti_spe)[1] <- "GENE"
```


```r
exp_conti_spe %>% write_tsv("MAGMA/conti_specificity_matrix_cluster83-205_excitatory-cortical-amygdala-hippocampal.txt")
```
