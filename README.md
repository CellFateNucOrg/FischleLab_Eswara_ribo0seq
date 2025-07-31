## RNAseq alignment and contrasts

Scripts for running nf_core pipelines in ./nf_core_scripts/ folder

Use ***./makeSampleSheet.R*** to create samples sheet and contrast file for running nf_core pipelines.

***nf_rnaseq_submit.sh*** runs the alignment pipeline on the cluster with STAR and salmon

***nf_differentialabundance.sh*** runs DESeq2 for specific contrasts

***gatherVolcanos.py*** pulls out volcano plots from nf_core subdirectories into a single pdf.


## Exploratory Analysis

Scripts in ./exploratory_analysis/ folder

***postProcessing.R*** annotates results tables with additional data from gtf (./custom/txt/\*_annotated.tsv) and
combines them into a single results file (./custom/rds/\*.results_annotated.RDS). Also produces a text file with number of
significant up/down regulated genes in each contras (./custom/txt/summaryUpDownDiffThresolds.txt)

***basicExploratoryPlots.R*** produces:
- static and interactive volcano plots at different zooms with more info than the pipeline volcanos (./custom/plots/volcano/)
- boxplots and barplots for LFC and up/down counts by chromosome (./custom/plots/byChromosome/)
- boxplots and barplots for LFC and up/down counts by chromosome region type (tip, arm, center) (./custom/plots/byChrRegion/)
- heatmaps of significantly changing genes (./custom/plots/heatmaps/)

***averageBigwigs.sh*** produces:
- average of bigwigs created by STAR for forward strand from all replicates 
(note: STAR labelling of forward and revese is backwards) (./custom/tracks/bigwigs_fr/\*.forward.avr.bw)
- average of bigwigs created by STAR for reverse strand from all replicates  (./custom/tracks/bigwigs_fr/\*.reverse.avr.bw)
- avereage of average forward and average reverse strand bigwigs to create
single track for each sample  (./custom/tracks/bigwigs_fr/\*.average_fr.bw)
- log2 fold change bigwigs of the different samples (./custom/tracks/bigwigs_compare/\*.average_fr.log2fc.bw)

***pairwiseComparisons.sh*** produces:
- correlation plots for pairs of samples
- venn diagrams for up/down regulated genes
- euler diagrams for triplets of up/down regulated genes

