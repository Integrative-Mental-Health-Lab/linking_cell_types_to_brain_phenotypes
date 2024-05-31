#!/usr/bin/bash
#SBATCH --time=2:00:00
#SBATCH --mem=16G

ml gcc
folder="SCZ"
file="PGC3_SCZ_wave3.european.autosome.public.v3.vcf.tsv.no_heading"
snp_col="ID"
p_col="PVAL"
ncol="NEFF"

MAGMA_v1.10/magma --annotate window=35,10 \
--snp-loc ${folder}/snploc_${file} \
--gene-loc auxfiles/NCBI37.3.gene.loc \
--out ${folder}/${file}.step1

MAGMA_v1.10/magma --bfile auxfiles/g1000_eur --pval ${folder}/${file} use=${snp_col},${p_col} ncol=${ncol} --gene-annot ${folder}/${file}.step1.genes.annot --out ${folder}/${file}.step2
