#!/bin/bash
#SBATCH --time=0-24:00:00
#SBATCH --mem-per-cpu=4G
#SBATCH --ntasks=1


source $CONDA_ACTIVATE env_nf

# percentages
export NXF_JVM_ARGS="-XX:InitialRAMPercentage=25 -XX:MaxRAMPercentage=75"

genomeVer=WS295
genomeDir=/mnt/external.data/MeisterLab/publicData/genomes/${genomeVer}
gtfFile=${genomeDir}/c_elegans.PRJNA13758.${genomeVer}.canonical_geneset_noRR_noSP.gtf

WORK_DIR=/mnt/external.data/MeisterLab/FischleLab_KarthikEswara/ribo0seq
CONFIG_FILE=/mnt/external.data/MeisterLab/nf-core/unibe_izb.config
#SRR_FILE=${WORK_DIR}/SRR_file.csv

nextflow run nf-core/differentialabundance -profile rnaseq,singularity --input ${WORK_DIR}/samplesheet.csv --outdir $WORK_DIR/diff_abund_2_canonical_noRRnoSP_moreSeq -c $CONFIG_FILE \
	--gtf $gtfFile   \
	--matrix ${WORK_DIR}/star_salmon/salmon.merged.gene_counts_noRR_noSP_moreSeq.tsv \
	--transcript_length_matrix ${WORK_DIR}/star_salmon/salmon.merged.gene_lengths.tsv \
	--contrasts ${WORK_DIR}/contrasts.csv \
	--deseq2_shrink_lfc true  --filtering_min_abundance 5 --filtering_min_samples 3
