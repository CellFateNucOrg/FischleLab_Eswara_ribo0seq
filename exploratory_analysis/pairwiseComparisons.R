library(GenomicRanges)
library(rtracklayer)
library(ggplot2)
library(dplyr)
library(BSgenome.Celegans.UCSC.ce11)
library(ggpubr)
library(plotly)
library(ggrepel)
library(rstatix)
library(tidyr)
library(gridExtra)
library(RColorBrewer)


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
