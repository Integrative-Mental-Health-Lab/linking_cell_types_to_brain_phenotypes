#!/usr/bin/bash
#SBATCH --time=2:00:00
#SBATCH --mem=16G

ml gcc
folder="MultipleSclerosis"
file="discovery_metav3.0.meta"
snp_col="SNP"
p_col="P"
N=41505

MAGMA_v1.10/magma --annotate window=35,10 \
--snp-loc ${folder}/snploc_${file} \
--gene-loc auxfiles/NCBI37.3.gene.loc \
--out ${folder}/${file}.step1

MAGMA_v1.10/magma --bfile auxfiles/g1000_eur --pval ${folder}/${file} use=${snp_col},${p_col} N=${N} --gene-annot ${folder}/${file}.step1.genes.annot --out ${folder}/${file}.step2
