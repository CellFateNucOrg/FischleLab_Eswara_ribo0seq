#!/bin/bash
#SBATCH --time=1-00:00:00
#SBATCH --mem-per-cpu=4G
#SBATCH --ntasks=1


source $CONDA_ACTIVATE env_nf

# percentages
export NXF_JVM_ARGS="-XX:InitialRAMPercentage=25 -XX:MaxRAMPercentage=75"

genomeVer=WS295
genomeDir=/mnt/external.data/MeisterLab/publicData/genomes/${genomeVer}
genomeFile=${genomeDir}/c_elegans.PRJNA13758.${genomeVer}.genomic.fa
gtfFile=${genomeDir}/c_elegans.PRJNA13758.${genomeVer}.canonical_geneset.gtf

WORK_DIR=/mnt/external.data/MeisterLab/FischleLab_KarthikEswara/ribo0seq
CONFIG_FILE=/mnt/external.data/MeisterLab/nf-core/unibe_izb.config


nextflow run nf-core/rnaseq -profile singularity -r 3.19.0 --input ${WORK_DIR}/samplesheet.csv --multiqc_title multiqc_rnaseq --outdir $WORK_DIR -c $CONFIG_FILE --fasta $genomeFile --gtf $gtfFile --pseudo_aligner salmon --save_unaligned -resume
