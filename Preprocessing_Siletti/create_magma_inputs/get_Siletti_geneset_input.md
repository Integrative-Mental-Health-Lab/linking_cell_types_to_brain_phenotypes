---
title: "Siletti 2022 scRNAseq Dataset - Get Gene Sets"
author: "Tayden Li"
output:
  html_document:
    keep_md: yes
---

Code is adapted from [Bryois (2020)](https://github.com/jbryois/scRNA_disease/tree/master).
The script produces MAGMA and LDSC inputs.

# Load Data

The Siletti data was downloaded from [here](https://github.com/linnarsson-lab/adult-human-brain).

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
file="/oak/stanford/groups/laramied/LDSC/Bryois2020/Code_Paper/Data/Siletti/adult_human_20221007.agg.loom"
h5f <- H5Fopen(file)
exp <- as.data.frame(t(h5f$matrix))
colnames(exp) <- h5f$col_attrs$Clusters
exp$Gene <- sub("\\..*", "", h5f$row_attrs$Accession)
h5closeAll()
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
exp_lvl5 <- transform(exp,ClusterID=paste0("Cluster", exp$column)) %>%
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
  mutate(Expr_sum_mean=Expr_sum_mean*1e6/sum(Expr_sum_mean))
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
n_genes_to_keep <- (n_genes * 0.10) %>% round()
```

### Functions

#### Get MAGMA input top10%


```r
magma_top10 <- function(d,Cell_type){
  d_spe <- d %>% group_by_(Cell_type) %>% top_n(.,n_genes_to_keep,specificity) 
  d_spe %>% do(write_group_magma(.,Cell_type))
}
```


```r
write_group_magma  = function(df,Cell_type) {
  df <- select(df,ClusterID,ENTREZ)
  df_name <- make.names(unique(df[1]))
  colnames(df)[2] <- df_name  
  dir.create(paste0("MAGMA/"), showWarnings = FALSE)
  select(df,2) %>% t() %>% as.data.frame() %>% rownames_to_column("Cat") %>%
  write_tsv("MAGMA/top10.txt",append=T)
return(df)
}
```

#### Get LDSC input top 10%


```r
write_group  = function(df,Cell_type) {
  df <- select(df,ClusterID,chr,start,end,ENTREZ)
  dir.create(paste0("LDSC/Bed"), showWarnings = FALSE,recursive = TRUE)
  write_tsv(df[-1],paste0("LDSC/Bed/",make.names(unique(df[1])),".bed"),col_names = F)
return(df)
}
```


```r
ldsc_bedfile <- function(d,Cell_type){
  d_spe <- d %>% group_by_(Cell_type) %>% top_n(.,n_genes_to_keep,specificity) 
  d_spe %>% do(write_group(.,Cell_type))
}
```

### Write MAGMA/LDSC input files 

Filter out genes with expression below 1 TPM.


```r
exp_lvl5 %>% filter(Expr_sum_mean>1) %>% magma_top10("ClusterID")
```

```
## Warning: `group_by_()` was deprecated in dplyr 0.7.0.
## ℹ Please use `group_by()` instead.
## ℹ See vignette('programming') for more help
## ℹ The deprecated feature was likely used in the dplyr package.
##   Please report the issue at <https://github.com/tidyverse/dplyr/issues>.
```

```
## # A tibble: 803,523 × 462
## # Groups:   ClusterID [461]
##    ClusterID Cluster0 Cluster1 Cluster10 Clust…¹ Clust…² Clust…³ Clust…⁴ Clust…⁵
##    <chr>     <chr>    <chr>    <chr>     <chr>   <chr>   <chr>   <chr>   <chr>  
##  1 Cluster0  146429   <NA>     <NA>      <NA>    <NA>    <NA>    <NA>    <NA>   
##  2 Cluster0  643866   <NA>     <NA>      <NA>    <NA>    <NA>    <NA>    <NA>   
##  3 Cluster0  9447     <NA>     <NA>      <NA>    <NA>    <NA>    <NA>    <NA>   
##  4 Cluster0  7161     <NA>     <NA>      <NA>    <NA>    <NA>    <NA>    <NA>   
##  5 Cluster0  8407     <NA>     <NA>      <NA>    <NA>    <NA>    <NA>    <NA>   
##  6 Cluster0  8556     <NA>     <NA>      <NA>    <NA>    <NA>    <NA>    <NA>   
##  7 Cluster0  55289    <NA>     <NA>      <NA>    <NA>    <NA>    <NA>    <NA>   
##  8 Cluster0  10451    <NA>     <NA>      <NA>    <NA>    <NA>    <NA>    <NA>   
##  9 Cluster0  131096   <NA>     <NA>      <NA>    <NA>    <NA>    <NA>    <NA>   
## 10 Cluster0  22861    <NA>     <NA>      <NA>    <NA>    <NA>    <NA>    <NA>   
## # … with 803,513 more rows, 453 more variables: Cluster105 <chr>,
## #   Cluster106 <chr>, Cluster107 <chr>, Cluster108 <chr>, Cluster109 <chr>,
## #   Cluster11 <chr>, Cluster110 <chr>, Cluster111 <chr>, Cluster112 <chr>,
## #   Cluster113 <chr>, Cluster114 <chr>, Cluster115 <chr>, Cluster116 <chr>,
## #   Cluster117 <chr>, Cluster118 <chr>, Cluster119 <chr>, Cluster12 <chr>,
## #   Cluster120 <chr>, Cluster121 <chr>, Cluster122 <chr>, Cluster123 <chr>,
## #   Cluster124 <chr>, Cluster125 <chr>, Cluster126 <chr>, Cluster127 <chr>, …
```


```r
exp_lvl5 %>% filter(Expr_sum_mean>1) %>% ldsc_bedfile("ClusterID")
```

```
## # A tibble: 803,523 × 5
## # Groups:   ClusterID [461]
##    ClusterID chr       start       end ENTREZ
##    <chr>     <chr>     <dbl>     <dbl> <chr> 
##  1 Cluster0  chr16  89162169  89367757 146429
##  2 Cluster0  chr14  24795738  24998731 643866
##  3 Cluster0  chr1  158928790 159146685 9447  
##  4 Cluster0  chr1    3469129   3752765 7161  
##  5 Cluster0  chr1  159787897 159995332 8407  
##  6 Cluster0  chr1  100710598 101085833 8556  
##  7 Cluster0  chr2  111390150 111975799 55289 
##  8 Cluster0  chr1  108013782 108607545 10451 
##  9 Cluster0  chr3   19089946  19677135 131096
## 10 Cluster0  chr17   5304719   5587832 22861 
## # … with 803,513 more rows
```
