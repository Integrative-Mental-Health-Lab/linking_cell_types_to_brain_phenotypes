#!/usr/bin/bash
#SBATCH --time=1:00:00
#SBATCH --mem=16G

ml gcc
folder="SCZ"
file="PGC3_SCZ_wave3.european.autosome.public.v3.vcf.tsv.no_heading"
gwas_name="SCZ_2022"
covar_file_spe="Siletti_l2_specificity_matrix_sig-only"
out_file_spe="Siletti_l2_joint_spe_sig-only"

MAGMA_v1.10/magma --gene-results ${folder}/${file}.step2.genes.raw --gene-covar ${folder}/${covar_file_spe}_${gwas_name}.txt --model joint-pairs direction=greater --out ${folder}/${out_file_spe}_${gwas_name}

