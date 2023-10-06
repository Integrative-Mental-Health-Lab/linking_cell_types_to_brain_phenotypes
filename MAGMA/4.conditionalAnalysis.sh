#!/usr/bin/bash
#SBATCH --time=1:00:00
#SBATCH --mem=16G

ml gcc
folder="MultipleSclerosis"
file="discovery_metav3.0.meta"
gwas_name="MultipleSclerosis_2018"
covar_file_spe="Siletti_l2_specificity_matrix_50-only"
out_file_spe="Siletti_l2_joint_spe_50-only"

MAGMA_v1.10/magma --gene-results ${folder}/${file}.step2.genes.raw --gene-covar ${folder}/${covar_file_spe}_${gwas_name}.txt --model joint-pairs direction=greater --out ${folder}/${out_file_spe}_${gwas_name}

