#!/usr/bin/bash
#SBATCH --time=1:00:00
#SBATCH --mem=16G

ml gcc
folder=MultipleSclerosis
file=discovery_metav3.0.meta
outfile_conti=Siletti_l2_conti-spe_MultipleSclerosis_2018
outfile_geneset=Siletti_l2_top10_MultipleSclerosis_2018

covar_file_conti="gene-level/Siletti_l2_conti_specificity_matrix.txt"
covar_file_geneset="gene-level/Siletti_l2_top10.txt"

MAGMA_v1.10/magma --gene-results ${folder}/${file}.step2.genes.raw --model direction=greater --gene-covar ${covar_file_conti} --out ${folder}/${outfile_conti}

MAGMA_v1.10/magma --gene-results ${folder}/${file}.step2.genes.raw --set-annot ${covar_file_geneset} --out ${folder}/${outfile_geneset}
