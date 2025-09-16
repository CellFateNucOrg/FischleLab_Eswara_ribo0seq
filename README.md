## RNAseq alignment and contrasts

Scripts for running nf_core pipelines in ./nf_core_scripts/ folder

Use ***./makeSampleSheet.R*** to create samples sheet and contrast file for running nf_core pipelines.

***nf_rnaseq_submit.sh*** runs the alignment pipeline on the cluster with STAR and salmon

***nf_differentialabundance.sh*** runs DESeq2 for specific contrasts

***gatherVolcanos.py*** pulls out volcano plots from nf_core subdirectories into a single pdf.


## Exploratory Analysis

Scripts in ./scripts_exploratory_analysis/ folder

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

***datasetComparisons.sh*** produces:
- correlation plots for pairs of samples (./custom/plots/correlation)
- upset plots for all N2 contrasts and lin-61 contrasts (./custom/plots/upset/)
- venn diagrams for up/down regulated genes for lin-61 contrasts (4-way) (./custom/plots/venn/)
- euler diagrams for triplets of up/down regulated genes fir kub-61 contrasts (all possible triplets) (./custom
)

***exploreChIP.R*** produces:
- bedfiles of singificantly up/down regulated genes from each mutant to be used by deeptools (./custom/exploreChIP/bed/)
- boxplots of LFC of genes from different datasets that overlap stratified 10kb triple peak regions, plotted either for all genes or only significantly upregulated autosomal arm genes (./custom/exploreChIP/*triplPeaks10kb_boxplot.pdf)
- histograms of the number of peaks per gene from each of the three peak-bedfiles plotted for all genes (./custom/exploreChIP/allGenes_bedPeakOverlap_histogram.pdf)
- stacked barplots showing the number of differnt types of peaks overlapped per gene for
triple H3K9me2/hpl2/lin61 combinations or just H3K9me2/hpl2 combinations, plotted for all genes or only significantly up/down regulated autosomal arm genes. (./custom/exploreChIP/bedPeakTypeOverlap_barplot.pdf)

***deeptools_sigGenes_v_ChIP.sh*** produces:
- deeptools heatmaps+profiles for ChIPseq of all histone marks over gene bodies +- 1kb 
for significantly up/down regulated genes in all four mutant contrasts. (./custom/exploreChIP/\*publicChIPseq_heatmap_sort\*.pdf)



## Final figures

Scripts in ./scripts_finalFigures/ folder

Focussing on significant genes (padj<0.05, LFC>0.5) on autosomal chromosome arms.

Using only 4 contrasts:

Only Dimerisation, Only LLPS, No LLPS No Dimerization and Double mutant compared to lin-61; hpl-2::gfp control.

***upregulatedGenesOnArms.R*** produces:
- Heatmap for all genes significantly upregulated in at least one group (combined or split by group)
- Boxplots of LFC of upregulated genes (or all expressed autosomal arm genes)
- Upset plot for all datasets
- Venn diagram for all datasets
- Pairwise Euler diagrams for all datasets
- Gene lists for significant genes in each group for GO and GSEA

output is all in ./custom/finalFigures/upregulatedOnArms/

***integrationWithChIP.R*** produces:
- euler plot for overlap between genes with peaks in different ChIP seq datasets (lfcByGroup_autosomalArmTripleChIPpeaks_padj0.05_lfc0.5.pdf)
- boxplot of LFC expression for core triple-peak genes in each mutant RNAseq (lfcByGroup_autosomalArmTripleChIPpeaks_padj0.05_lfc0.5.pdf)
- same as above but with ChIP seq signal stratified into three quantiles (Q1 lowwest-Q3 highest) (lfcByGroup_autosomalArmTripleChIPpeaks_padj0.05_lfc0.5_quantile.pdf)

output is all in ./custom/finalFigures/integrationWithChIP/
