# Plotting significant gene number and LFC by chromosome, chromosome region (arm/center),
# volcano plots and heatmaps

library(GenomicRanges)
library(rtracklayer)
library(ggplot2)
library(dplyr)
library(BSgenome.Celegans.UCSC.ce11)
library(ggpubr)
library(plotly)
library(ggrepel)
library(rstatix)
library(htmlwidgets)
library(RColorBrewer)
library(ComplexHeatmap)



theme_set(
  theme_classic()+
    theme(panel.grid.major = element_blank(),
          panel.grid.minor = element_blank(),
          axis.title.y=ggtext::element_markdown(size=9),
          axis.title.x=ggtext::element_markdown(size=9),
          title=ggtext::element_markdown(size=9)
    )
)

serverPath="/Volumes/external.data/MeisterLab"
#serverPath="Z:/MeisterLab"

workDir=paste0(serverPath,"/FischleLab_KarthikEswara/ribo0seq")
runName="/diff_abund_3_canonical_noRRnoSP_moreSeq"

source(paste0(workDir,"/scripts_finalFigures/functions_finalFigures.R"))

contrasts<-read.csv(paste0(workDir,"/contrasts.csv"),sep=",",header=T)

contrastsToKeep<-10:13
contrasts<-contrasts[contrastsToKeep,]
contrasts$prettyName<-c("Only Dimerization","Only LLPS", "No LLPS No Dimerization", "Double mutant")
write.csv(contrasts,paste0(workDir,"/contrasts_subset.csv"),row.names=F, quote=T)

contrasts<-read.csv(paste0(workDir,"/contrasts_subset.csv"),sep=",",header=T)

prefix="ribo0_canonical_geneset_all"

setwd(workDir)

dir.create(paste0(workDir,runName,"/custom/finalFigures"), showWarnings = FALSE, recursive = TRUE)



## subset data -------
results<-readRDS(paste0(workDir,runName,"/custom/rds/",prefix,".results_annotated.RDS"))
dir.create(paste0(workDir,runName,"/custom/finalFigures/upregulatedOnArms"), showWarnings = FALSE, recursive = TRUE)
rockman<-readRDS(paste0(serverPath,"/publicData/Various/chrRegions_Rockman2009.RDS"))
gr<-tableToGranges(results,sort=FALSE)
seqlevelsStyle(gr)<-"UCSC"
ol<-findOverlaps(resize(gr,width=1,fix="start"),rockman,ignore.strand=T)
gr$chrRegionType[queryHits(ol)] <- paste0(rockman$chrRegionType[subjectHits(ol)])
gr$chrRegion[queryHits(ol)] <- paste0(seqnames(rockman)[subjectHits(ol)],"_",rockman$chrRegionType[subjectHits(ol)])
res<-as.data.frame(gr)
contrasts

lfcVal=0.5
padjVal=0.05

ss<-res %>% filter(shortId %in% contrasts$shortId,
               chrRegionType %in% c("arm", "tip"),
               seqnames %in% seqnames(Celegans)[1:5],
               padj<padjVal,
               log2FoldChange>lfcVal)
ss$longId<-droplevels(ss$group)
ss$shortId<-droplevels(ss$shortId)
ss$group<-factor(ss$shortId,levels=contrasts$shortId,labels=contrasts$prettyName)
table(ss$seqnames,ss$group)
saveRDS(ss, paste0(workDir,runName,"/custom/rds/",prefix,"_subset.results_annotated.RDS"))

# keep values in all samples for genes significant in at least one
genesToKeep<-unique(ss$gene_id)
ss<-res[res$gene_id %in% genesToKeep & res$shortId %in% contrasts$shortId,]
ss$longId<-droplevels(ss$group)
ss$shortId<-droplevels(ss$shortId)
ss$group<-ss$shortId
ss$group<-factor(ss$shortId,levels=contrasts$shortId,labels=contrasts$prettyName)
table(ss$seqnames,ss$group)
saveRDS(ss, paste0(workDir,runName,"/custom/rds/",prefix,"_subsetByGene.results_annotated.RDS"))

rm(results)
rm(ss)
rm(genesToKeep)

## Heatmap of significant genes -----
res<-readRDS(paste0(workDir,runName,"/custom/rds/",prefix,"_subsetByGene.results_annotated.RDS"))

#makeHeatmapPlot(res, numSamplesSignificant=1, lfcVal, padjVal, direction="up",
#                groupsToPlot=levels(res$group), setName="mutants",
#                chromosomes="autosomal")

chromosomes<-"autosomal"
direction<-"up"
numSamplesSignificant=1

#mat_lfc<-getLFCMat(res, numSamplesSignificant, lfcVal, padjVal, direction, chromosomes)

mat_padj<-gatherResults(res,valueColumn="padj")
mat_padj[is.na(mat_padj)]<-1
mat_lfc<-gatherResults(res,valueColumn="log2FoldChange")
sig<-(mat_padj[,2:ncol(mat_padj)]<padjVal & mat_lfc[,2:ncol(mat_lfc)]>lfcVal)
sigGenes<-mat_padj$gene_id[rowSums(sig)>=numSamplesSignificant]
mat_lfc<-mat_lfc[mat_lfc$gene_id %in% sigGenes,]
rownames(mat_lfc)<-mat_lfc$gene_id
mat_lfc<-mat_lfc[,2:ncol(mat_lfc)]
mat_lfc<-as.matrix(mat_lfc)
colnames(mat_lfc)<-gsub("log2FoldChange_","",colnames(mat_lfc))

clusters<-data.frame(name=apply(sig,1,function(x) paste0(gsub("padj_","",names(x)[x]),collapse=" & ")),
                     rank=rowSums(sig))

groups<-clusters[!duplicated(clusters),]
groups<-groups[order(-groups$rank),]

clusters$name<-factor(clusters$name,levels=groups$name)

mat_list <- lapply(split(seq_len(nrow(mat_lfc)), clusters$name), function(i) mat_lfc[i, , drop = FALSE])

plotTitle<-paste0(nrow(mat_lfc)," ", chromosomes," arm genes sig. in \n>=",
                  numSamplesSignificant," sample (padj<",padjVal," ",
                  ifelse(direction=="both","|LFC|","LFC"),
                  ifelse(direction=="down","< -",">"),lfcVal,")")
# Shared color mapping
col_fun <- circlize::colorRamp2(breaks = c(-2, 0, 2), colors = c("blue", "white", "red"))

pdf(paste0(workDir, runName,"/custom/finalFigures/upregulatedOnArms/hclust_heatmap_sigSamples",
           numSamplesSignificant,"_padj",padjVal,"_lfc",lfcVal,"_",direction,"_",chromosomes,"Chr.pdf"),
    width=4,height=9)

ht<-Heatmap(mat_lfc,show_row_names=F,
            cluster_columns=F, cluster_rows=T,  show_row_dend = F,
            column_names_rot=90, column_title=plotTitle, row_title_rot=0,
            column_names_gp=gpar(fontsize=10),
            col=col_fun,
            column_title_gp = gpar(fontsize = 12, fontface = "bold", just = 0),
            heatmap_legend_param=list(title=expression("log"[2]~"(FC)")),
            )
ht<-draw(ht)
dev.off()


pdf(paste0(workDir, runName,"/custom/finalFigures/upregulatedOnArms/hclust_heatmap_sigSamples",
           numSamplesSignificant,"_padj",padjVal,"_lfc",lfcVal,"_",direction,"_",chromosomes,"Chr_grouped.pdf"),
    width=8,height=11)

# Create a heatmap for each matrix
ht_list <- lapply(names(mat_list), function(name) {
  Heatmap(mat_list[[name]],row_title=name, show_row_names=F,
          cluster_columns=F, cluster_rows=T,  show_row_dend = F,
          column_names_rot=45, column_title=plotTitle,
          column_names_gp = gpar(fontsize = 10),
          column_title_gp = gpar(fontsize = 10, fontface = "bold", just = -5),
          row_title_rot=0, row_title_gp=gpar(fontsize=10),
          col=col_fun,
          heatmap_legend_param=list(title=expression("log"[2]~"(FC)"))
  )
})

# Combine them into one plot
ht_combined <- Reduce(`%v%`, ht_list)
draw(ht_combined,merge_legend=T, heatmap_legend_side="right")
dev.off()

sapply(mat_list,nrow)
