library(GenomicRanges)
library(rtracklayer)
library(ggplot2)
library(dplyr)
library(BSgenome.Celegans.UCSC.ce11)
library(ggpubr)
#library(plotly)
#library(ggrepel)
#library(rstatix)
library(tidyr)
library(grid)
library(gridExtra)
library(ggplotify)
library(RColorBrewer)
library(eulerr)
library(ComplexUpset)
library(ggVennDiagram)


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
runName="/diff_abund_2_canonical_noRRnoSP"

contrasts<-read.csv(paste0(workDir,"/contrasts.csv"),sep=",",header=T)
prefix="ribo0_canonical_geneset_all"

setwd(workDir)

colorSet<-brewer.pal(8, "Dark2")


## functions ------
gatherResults<-function(results, valueColumn, nameColumn="group"){
  results[,nameColumn]<-droplevels(results[,nameColumn])
  for(c in levels(results[,nameColumn])){
    tmp<-results[results[,nameColumn]==c,c("gene_id",valueColumn)]
    colnames(tmp)[2]<-paste0(valueColumn,"_",c)
    if(!exists("gathered")){
      gathered<-tmp
    } else {
      gathered<-merge(gathered,tmp,by="gene_id",all=T)
    }
  }
  return(gathered)
}

#' Get a binary matrix of significant results by group
#'
#' From a results table, extract genes that are signficant in at least
#' `numSamplesSignificant` samples, with an absolute log2 fold change greater
#' than `lfcVal` and adjusted p-value less than `padjVal`. Log2FoldChange values
#' can be evaluated in both directions (up and down), or just up or down
#' regulated.
#' @param results A data frame containing the results of differential expression analysis.
#' @param numSamplesSignificant Minimum number of samples in which a gene must be significant.
#' @param lfcVal Minimum absolute log2 fold change value for significance.
#' @param padjVal Maximum adjusted p-value for significance.
#' @param direction Direction of regulation to consider: "both", "up", or "down".
#' @param chromosomes Chromosomes to be considered in the analysis: "all" or "autosomal".
getBinaryMat<-function(results, numSamplesSignificant=1, lfcVal=0, padjVal=0.05,
                       direction="both", chromosomes="all"){
  if(chromosomes=="all"){
    chr<-c("I","II","III","IV","V","X")
    res<-results[results$seqnames %in% chr,]
  } else if(chromosomes=="autosomal"){
    chr<-c("I","II","III","IV","V")
    res<-results[results$seqnames %in% chr,]
  } else {
    stop("Invalid chromosomes specified. Use 'all' or 'autosomal'.")
  }
  mat_padj<-gatherResults(res,valueColumn="padj")
  mat_padj[is.na(mat_padj)]<-1
  mat_lfc<-gatherResults(res,valueColumn="log2FoldChange")
  if(direction=="both"){
    sig<-(mat_padj[,2:ncol(mat_padj)]<padjVal & abs(mat_lfc[,2:ncol(mat_lfc)])>lfcVal)
  } else if(direction=="up"){
    sig<-(mat_padj[,2:ncol(mat_padj)]<padjVal & mat_lfc[,2:ncol(mat_lfc)]>lfcVal)
  } else if(direction=="down"){
    sig<-(mat_padj[,2:ncol(mat_padj)]<padjVal & mat_lfc[,2:ncol(mat_lfc)]<(-lfcVal))
  } else {
    stop("Invalid direction specified. Use 'both', 'up', or 'down'.")
  }
  rownames(sig)<-mat_padj[,1]
  colnames(sig)<-gsub("padj_","",colnames(sig))
  sig<-sig[rowSums(sig)>=numSamplesSignificant,]

  return(sig)
}


## correlations -----
# N2 contrasts
results<-readRDS(paste0(workDir, runName,"/custom/rds/",prefix,".results_annotated.RDS"))
dir.create(paste0(workDir, runName,"/custom/plots/correlation"), showWarnings = FALSE, recursive = TRUE)

n2contrasts<-grep("_vs_N2$",levels(results$group),value=T)
lin61contrasts<-grep("_vs_HPL2GFP__lin61$",levels(results$group),value=T)

n2pairs<-combn(n2contrasts,2)
lin61pairs<-combn(lin61contrasts,2)
n2pairs
lin61pairs

group_colors <- setNames(colorSet, n2contrasts)

#c=1
plotList=list()
for(c in 1:ncol(n2pairs)){
  c1<-n2pairs[1,c]
  c2<-n2pairs[2,c]

  df<-results %>%
    dplyr::filter(group %in% c(c1,c2)) %>%
    dplyr::select(gene_id, group, log2FoldChange) %>%
    tidyr::pivot_wider(names_from=group, values_from=c(log2FoldChange))

  corVal<-cor(df[,c1], df[,c2])

  p<-ggplot(df, aes(x=.data[[c1]], y=.data[[c2]])) +
    geom_point(size=0.5,alpha=0.3) +
    geom_smooth(method="lm", color="blue", se=FALSE) +
    coord_cartesian(xlim=c(-4,4), ylim=c(-4,4)) +
    geom_hline(yintercept=0, linetype="dashed", color="grey") +
    geom_vline(xintercept=0, linetype="dashed", color="grey") +
    xlab(c1) + ylab(c2)+
    ggpubr::stat_cor(aes(label = ..r.label..), method="pearson",
                     cor.coef.name = c("R"), output.type = "text",
                     label.x=-4, label.y=4) +
    ggpubr::stat_cor(aes(label = ..r.label..), method="spearman",
                   cor.coef.name = "\U03C1", output.type = "text",
                   label.x=-4, label.y=3.5)+
    theme(axis.title.x = ggtext::element_markdown(colour = as.character(group_colors[c1])),
          axis.title.y = ggtext::element_markdown(colour = as.character(group_colors[c2])))
  p
  plotList[[paste0(c1,"__",c2)]]<-p
}

rowNum=3
colNum=4
if(length(plotList) < rowNum*colNum){
  rowNum=ceiling(length(plotList)/colNum)
}
numPages<-ceiling(length(plotList)/rowNum/colNum)

for(i in 1:numPages){
  startIndex<-(i-1)*rowNum*colNum+1
  endIndex<-min(i*rowNum*colNum, length(plotList))
  plotListPage<-plotList[startIndex:endIndex]

  p<-ggarrange(plotlist=plotListPage, ncol=colNum, nrow=rowNum, common.legend = TRUE, legend="bottom")

  ggsave(paste0(workDir, runName,"/custom/plots/correlation/n2_contrasts_correlation_page",i,".png"), plot=p, width=12, height=8)
}



group_colors <- setNames(colorSet[1:length(lin61contrasts)], lin61contrasts)

plotList=list()
for(c in 1:ncol(lin61pairs)){
  c1<-lin61pairs[1,c]
  c2<-lin61pairs[2,c]

  df<-results %>%
    dplyr::filter(group %in% c(c1,c2)) %>%
    dplyr::select(gene_id, group, log2FoldChange) %>%
    tidyr::pivot_wider(names_from=group, values_from=c(log2FoldChange))

  corVal<-cor(df[,c1], df[,c2])

  p<-ggplot(df, aes(x=.data[[c1]], y=.data[[c2]])) +
    geom_point(size=0.5,alpha=0.3) +
    geom_smooth(method="lm", color="blue", se=FALSE) +
    coord_cartesian(xlim=c(-4,4), ylim=c(-4,4)) +
    geom_hline(yintercept=0, linetype="dashed", color="grey") +
    geom_vline(xintercept=0, linetype="dashed", color="grey") +
    xlab(c1) + ylab(c2)+
    ggpubr::stat_cor(aes(label = ..r.label..), method="pearson",
                     cor.coef.name = c("R"), output.type = "text",
                     label.x=-4, label.y=4) +
    ggpubr::stat_cor(aes(label = ..r.label..), method="spearman",
                     cor.coef.name = "\U03C1", output.type = "text",
                     label.x=-4, label.y=3.5)+
    theme(axis.title.x = ggtext::element_markdown(colour = as.character(group_colors[c1])),
          axis.title.y = ggtext::element_markdown(colour = as.character(group_colors[c2])))
  p
  plotList[[paste0(c1,"__",c2)]]<-p
}

rowNum=3
colNum=4
if(length(plotList) < rowNum*colNum){
  rowNum=ceiling(length(plotList)/colNum)
}
numPages<-ceiling(length(plotList)/rowNum/colNum)

for(i in 1:numPages){
  startIndex<-(i-1)*rowNum*colNum+1
  endIndex<-min(i*rowNum*colNum, length(plotList))
  plotListPage<-plotList[startIndex:endIndex]

  p<-ggarrange(plotlist=plotListPage, ncol=colNum, nrow=rowNum, common.legend = TRUE, legend="bottom")

  ggsave(paste0(workDir, runName,"/custom/plots/correlation/lin61_contrasts_correlation_page",i,".png"), plot=p, width=12, height=8)
}


## upset plots -------
results<-readRDS(paste0(workDir,runName,"/custom/rds/",prefix,".results_annotated.RDS"))
dir.create(paste0(workDir, runName, "/custom/plots/upset"), showWarnings = FALSE, recursive = TRUE)
n2contrasts<-grep("vs_N2$",levels(results$group),value=T)
lin61contrasts<-grep("vs_HPL2GFP__lin61$",levels(results$group),value=T)

res_n2<-results[results$group %in% n2contrasts,]
res_lin61<-results[results$group %in% lin61contrasts[2:5],]

#' Make upset plot
#'
#' Make an upset plot for genes in `res` that are significant in at least `numSamplesSignificant`
#' groups, using a log2 fold change threshold of `lfcVal` and an adjusted p-value threshold of `padjVal`.
#' @param res Results data frame containing gene expression data.
#' @param numSamplesSignificant Minimum number of samples in which a gene must be significant.
#' @param lfcVal Minimum absolute log2 fold change value for significance.
#' @param padjVal Maximum adjusted p-value for significance.
#' @param direction Direction of regulation to consider: "both", "up", or "down".
#' @param minSize Minimum size of the sets to be included in the upset plot.
#' @param minDegree Minimum degree of the sets to be included in the upset plot.
#' @param groupsToPlot Names of the groups to be included in the upset plot. If NULL, all groups are included.
#' @param setName Name of the set to be included in the plot filename.
#' @param chromosomes Chromosomes to be considered in the analysis: "all" or "autosomal".
#' @return A ggplot object representing the upset plot. it is also saved as a PNG file in the ./custom/plots/upset/.
makeUpsetPlot<-function(res, numSamplesSignificant, lfcVal, padjVal, direction,
                        minSize, minDegree, groupsToPlot=NULL, setName="",
                        chromosomes="all"){
  binMat<-getBinaryMat(res, numSamplesSignificant, lfcVal, padjVal, direction, chromosomes)
  if(is.null(groupsToPlot)){
    groupsToPlot<-colnames(binMat)
  }
  if(!(setName=="")){
    setName=paste0(setName,"_")
  }
  p<-ComplexUpset::upset(data.frame(binMat[,groupsToPlot]),groupsToPlot,
                         min_size=minSize, min_degree=minDegree,
                         base_annotations=list(
                           'Intersection ratio'=intersection_ratio(text_mapping=aes(label=!!upset_text_percentage())),
                           'Intersection size'=intersection_size()))
  ggsave(paste0(workDir, runName,"/custom/plots/upset/",setName,"upsetPlot_minSize",minSize,"_minDeg",minDegree,"_",direction,"_",chromosomes,"Chr_lfcVal",lfcVal,".png"), plot=p, width=12, height=8)
  return(p)
}


numSamplesSignificant=1
padjVal=0.05


## N2  both directions
makeUpsetPlot(res_n2, numSamplesSignificant, lfcVal=1, padjVal, direction="both",
                        minSize=50, minDegree=2, groupsToPlot=n2contrasts, setName="n2_contrasts",
                        chromosomes="all")

makeUpsetPlot(res_n2, numSamplesSignificant, lfcVal=1, padjVal, direction="both",
              minSize=20, minDegree=3, groupsToPlot=n2contrasts, setName="n2_contrasts",
              chromosomes="all")

# N2 upregulated
makeUpsetPlot(res_n2, numSamplesSignificant, lfcVal=1, padjVal, direction="up",
              minSize=20, minDegree=2, groupsToPlot=n2contrasts, setName="n2_contrasts",
              chromosomes="all")

# N2 downregulated
makeUpsetPlot(res_n2, numSamplesSignificant, lfcVal=1, padjVal, direction="down",
              minSize=20, minDegree=2, groupsToPlot=n2contrasts, setName="n2_contrasts",
              chromosomes="all")

# N2 upregulated autosomes
makeUpsetPlot(res_n2, numSamplesSignificant, lfcVal=1, padjVal, direction="up",
              minSize=20, minDegree=2, groupsToPlot=n2contrasts, setName="n2_contrasts",
              chromosomes="autosomal")

# lower LFCval
makeUpsetPlot(res_n2, numSamplesSignificant, lfcVal=0.5, padjVal, direction="up",
              minSize=20, minDegree=2, groupsToPlot=n2contrasts, setName="n2_contrasts",
              chromosomes="autosomal")

makeUpsetPlot(res_n2, numSamplesSignificant, lfcVal=0, padjVal, direction="up",
              minSize=20, minDegree=2, groupsToPlot=n2contrasts, setName="n2_contrasts",
              chromosomes="autosomal")



minSize=1

## lin-61  both directions
makeUpsetPlot(res_lin61, numSamplesSignificant, lfcVal=1, padjVal, direction="both",
              minSize=minSize, minDegree=2, groupsToPlot=lin61contrasts[2:5], setName="lin61_contrasts",
              chromosomes="all")

# lin-61 upregulated
makeUpsetPlot(res_lin61, numSamplesSignificant, lfcVal=1, padjVal, direction="up",
              minSize=minSize, minDegree=2, groupsToPlot=lin61contrasts[2:5], setName="lin61_contrasts",
              chromosomes="all")

# lin-61 downregulated
makeUpsetPlot(res_lin61, numSamplesSignificant, lfcVal=1, padjVal, direction="down",
              minSize=minSize, minDegree=2, groupsToPlot=lin61contrasts[2:5], setName="lin61_contrasts",
              chromosomes="all")



#  lin-61 upregulated autosomes
makeUpsetPlot(res_lin61, numSamplesSignificant, lfcVal=1, padjVal, direction="up",
              minSize=minSize, minDegree=2, groupsToPlot=lin61contrasts[2:5], setName="lin61_contrasts",
              chromosomes="autosomal")

# lin-61 lower LFCval
makeUpsetPlot(res_lin61, numSamplesSignificant, lfcVal=0.5, padjVal, direction="up",
              minSize=minSize, minDegree=2, groupsToPlot=lin61contrasts[2:5], setName="lin61_contrasts",
              chromosomes="autosomal")

makeUpsetPlot(res_lin61, numSamplesSignificant, lfcVal=0, padjVal, direction="up",
              minSize=minSize, minDegree=2, groupsToPlot=lin61contrasts[2:5], setName="lin61_contrasts",
              chromosomes="autosomal")



## euler plots ------
results<-readRDS(paste0(workDir,runName,"/custom/rds/",prefix,".results_annotated.RDS"))
dir.create(paste0(workDir, runName, "/custom/plots/euler"), showWarnings = FALSE, recursive = TRUE)
dir.create(paste0(workDir, runName, "/custom/plots/venn"), showWarnings = FALSE, recursive = TRUE)

lin61contrasts<-grep("vs_HPL2GFP__lin61$",levels(results$group),value=T)
res_lin61<-results[results$group %in% lin61contrasts[2:5],]


#' Make venn and euler plots f
#'
#' Make an venn/euler plots for genes in `res` that are significant in at least `numSamplesSignificant`
#' groups, using a log2 fold change threshold of `lfcVal` and an adjusted p-value threshold of `padjVal`.
#' @param res Results data frame containing gene expression data.
#' @param numSamplesSignificant Minimum number of samples in which a gene must be significant.
#' @param lfcVal Minimum absolute log2 fold change value for significance.
#' @param padjVal Maximum adjusted p-value for significance.
#' @param direction Direction of regulation to consider: "both", "up", or "down".
#' @param groupsToPlot Names of the groups to be included in the upset plot. If NULL, all groups are included.
#' @param setName Name of the set to be included in the plot filename.
#' @param chromosomes Chromosomes to be considered in the analysis: "all" or "autosomal".
#' @return A ggplot object representing the upset plot. it is also saved as a PNG file in the ./custom/plots/upset/.
makeEulerVennPlots<-function(res, numSamplesSignificant=1, lfcVal, padjVal, direction,
                        groupsToPlot=NULL, setName="",chromosomes="all"){
  binMat<-getBinaryMat(res, numSamplesSignificant, lfcVal, padjVal, direction, chromosomes)
  if(is.null(groupsToPlot)){
    groupsToPlot<-colnames(binMat)
  }
  if(!(setName=="")){
    setName=paste0(setName,"_")
  }
  ### 4-way
  geneList<-list()
  for(c in lin61contrasts[2:5]){
    geneList[[c]]<-rownames(binMat[binMat[,c]==1,])
  }

  p1<-ggVennDiagram(geneList, label_alpha = 0) +
    scale_fill_gradient(low = "white", high="darkred") +
    coord_cartesian(clip="off") +
    theme(plot.margin = margin(1, 1, 1, 3, "cm"))
  ggsave(paste0(workDir, runName, "/custom/plots/venn/",setName,"_4way_combination_",direction,"_",chromosomes,"Chr_lfcVal",lfcVal,".pdf"), plot=p1, width=8, height=6)

  ### 3-way combinations
  lin61groups<-combn(lin61contrasts[2:5],3)
  lin61groups
  colors <- brewer.pal(4, "Accent")
  named_colors <- setNames(colors,lin61contrasts[2:5])
  for(i in 1:ncol(lin61groups)){
    fit<-euler(binMat[,lin61groups[,i]])
    #fit
    p<-plot(fit,quantities=T, fills = list(fill = named_colors[lin61groups[,i]], alpha = 0.7))
    ggsave(paste0(workDir, runName, "/custom/plots/euler/",setName,"_3way_",direction,"_",chromosomes,"Chr_lfcVal",lfcVal,"_combination",i,".pdf"), plot=p, width=8, height=6)
  }
  return(p1)
}


numSamplesSignificant=1
padjVal=0.05
#binMat<-getBinaryMat(res_lin61, numSamplesSignificant, lfcVal=1, padjVal, direction="both", chromosomes="all")


## lin61
makeEulerVennPlots(res_lin61, numSamplesSignificant, lfcVal=1, padjVal, direction="both",
                             groupsToPlot=lin61contrasts[2:5], setName="lin61_contrasts",
                             chromosomes="all")

makeEulerVennPlots(res_lin61, numSamplesSignificant, lfcVal=1, padjVal, direction="up",
                   groupsToPlot=lin61contrasts[2:5], setName="lin61_contrasts",
                   chromosomes="all")


makeEulerVennPlots(res_lin61, numSamplesSignificant, lfcVal=1, padjVal, direction="down",
                   groupsToPlot=lin61contrasts[2:5], setName="lin61_contrasts",
                   chromosomes="all")

makeEulerVennPlots(res_lin61, numSamplesSignificant, lfcVal=1, padjVal, direction="up",
                   groupsToPlot=lin61contrasts[2:5], setName="lin61_contrasts",
                   chromosomes="autosomal")

makeEulerVennPlots(res_lin61, numSamplesSignificant, lfcVal=0.5, padjVal, direction="up",
                   groupsToPlot=lin61contrasts[2:5], setName="lin61_contrasts",
                   chromosomes="autosomal")

makeEulerVennPlots(res_lin61, numSamplesSignificant, lfcVal=0, padjVal, direction="up",
                   groupsToPlot=lin61contrasts[2:5], setName="lin61_contrasts",
                   chromosomes="autosomal")

