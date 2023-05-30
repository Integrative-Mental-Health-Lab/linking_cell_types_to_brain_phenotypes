#!/usr/bin/bash
#SBATCH --time=1:00:00
#SBATCH --mem=16G

ml gcc
folder="PTSD"
file="eur_ptsd_pcs_v4_aug3_2021.allchrX.fuma"
gwas_name="PTSD_wave3"
covar_file_spe="Siletti_l2_specificity_matrix_50-only"
covar_file_geneset="Siletti_l2_top10_geneset_50-only"
out_file_spe="Siletti_l2_joint_spe_50-only"
out_file_geneset="Siletti_l2_joint_top10_50-only"

MAGMA_v1.10/magma --gene-results ${folder}/${file}.step2.genes.raw --gene-covar ${folder}/${covar_file_spe}_${gwas_name}.txt --model joint-pairs direction=greater --out ${folder}/${out_file_spe}_${gwas_name}

MAGMA_v1.10/magma --gene-results ${folder}/${file}.step2.genes.raw --set-annot ${folder}/${covar_file_geneset}_${gwas_name}.txt --model joint-pairs direction=greater --out ${folder}/${out_file_geneset}_${gwas_name}
