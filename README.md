## RNAseq alignment and contrasts
Running nf_core pipelines.
Use _./makeSampleSheet.R_ to create samples sheet and contrast file for running nf_core pipelines.
_nf_rnaseq_submit.sh_ runs the alignment pipeline on the cluster with STAR and salmon
_nf_differentialabundance.sh_ runs DESeq2 for specific contrasts

## Basic plots
_gatherVolcanos.py_ pulls out volcano plots from subdirectories into a single pdf.
_postProcessing.R_ produces: 
- static and interactive volcano plots at different zooms with more info than the pipeline volcanos
- boxplots and barplots for LFC and up/down counts by chromosome
- boxplots and barplots for LFC and up/down counts by chromosome region type (tip, arm, center)
