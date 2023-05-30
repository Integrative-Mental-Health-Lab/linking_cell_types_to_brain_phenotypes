#!/usr/bin/bash
#SBATCH --time=1:00:00
#SBATCH --mem=16G

ml gcc
folder=SCZ
file=PGC3_SCZ_wave3.european.autosome.public.v3.vcf.tsv.no_heading

outfile_conti=Siletti_l2_conti-spe_SCZ_2022
covar_file_conti="gene-level/Siletti_l2_conti_specificity_matrix.txt"
MAGMA_v1.10/magma --gene-results ${folder}/${file}.step2.genes.raw --model direction=greater --gene-covar ${covar_file_conti} --out ${folder}/${outfile_conti}

outfile_conti=Siletti_l2_conti-spe_cluster239-296_SCZ_2022
covar_file_conti="gene-level/Siletti_l2_spe_matrix_cluster239-296_interneurons.txt"
MAGMA_v1.10/magma --gene-results ${folder}/${file}.step2.genes.raw --model direction=greater --gene-covar ${covar_file_conti} --out ${folder}/${outfile_conti}


outfile_conti=Siletti_l2_conti-spe_cluster83-205_SCZ_2022
covar_file_conti="gene-level/Siletti_l2_spe_matrix_cluster83-205_excitatory-cortical-amygdala-hippocampal.txt"
MAGMA_v1.10/magma --gene-results ${folder}/${file}.step2.genes.raw --model direction=greater --gene-covar ${covar_file_conti} --out ${folder}/${outfile_conti}


outfile_conti=Siletti_l2_conti-spe_cluster83-460_SCZ_2022
covar_file_conti="gene-level/Siletti_l2_spe_matrix_cluster83-460_neurons.txt"
MAGMA_v1.10/magma --gene-results ${folder}/${file}.step2.genes.raw --model direction=greater --gene-covar ${covar_file_conti} --out ${folder}/${outfile_conti}

