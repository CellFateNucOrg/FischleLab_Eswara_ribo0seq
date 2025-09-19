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
runName="/diff_abund_3_canonical_noRRnoSP_moreSeq"

source(paste0(workDir,"/scripts_finalFigures/functions_finalFigures.R"))

contrasts<-read.csv(paste0(workDir,"/contrasts_subset.csv"),sep=",",header=T)

prefix="ribo0_canonical_geneset_all"

setwd(workDir)

dir.create(paste0(workDir,runName,"/custom/finalFigures/integrationWithChIP"), showWarnings = FALSE, recursive = TRUE)

lfcVal=0.5
padjVal=0.05

# peaks
bedFiles<-list.files(paste0(workDir,"/../GSE271919_ChIPseq/bed/"),pattern=".*GSE.*\\.bed")
targets<-sapply(strsplit(bedFiles,"_"),"[[",1)
chosen<-c(7,9,5)
beddf<-data.frame(bedFiles=paste0(workDir,"/../GSE271919_ChIPseq/bed/",bedFiles[chosen]),targets=targets[chosen])

# signal
bwFiles<-list.files(paste0(workDir,"/../GSE271919_ChIPseq/bigwig/"),pattern=".*input_normalised\\.bw")
targets<-sapply(strsplit(bwFiles,"_"),"[[",1)
chosen<-c(7,9,5)
bwdf<-data.frame(bwFiles=paste0(workDir,"/../GSE271919_ChIPseq/bigwig/",bwFiles[chosen]),targets=targets[chosen])


# genes
results<-readRDS(paste0(workDir,runName,"/custom/rds/",prefix,"_regionType.results_annotated.RDS"))
results$padj[is.na(results$padj)]<-1
results$longId<-droplevels(results$group)
results$shortId<-droplevels(results$shortId)
results$group<-factor(results$shortId,levels=contrasts$shortId,labels=contrasts$prettyName)
results<-results[results$group %in% contrasts$prettyName,]
results<-results[results$seqnames!="chrM",]

gr<-tableToGranges(results)

for(i in 1:nrow(beddf)){
  bed<-import(beddf$bedFiles[i])
  seqlevelsStyle(bed)<-"UCSC"
  cname<-paste0(beddf$target[i],"_peaks")
  mcols(gr)[cname]<-countOverlaps(gr,bed)
}

gr
peakCols<-grep("_peaks",colnames(mcols(gr)))
gr$peakTypeNumber<-rowSums(data.frame(mcols(gr)[,peakCols])>0)

bwdf$totalScore<-0
for(i in 1:nrow(bwdf)){
  bw<-import(bwdf$bwFiles[i])
  bwdf$totalScore[i]<-sum(bw$score)
  seqlevelsStyle(bw)<-"UCSC"
  seqlevels(bw)<-seqlevels(gr)
  cov<-coverage(bw,weight="score")
  cname<-paste0(bwdf$target[i],"_signal")
  print(cname)
  gr<-binnedAverage(gr,cov,varname=cname,na.rm=T)
}
i=1
# for(i in 1:nrow(bwdf)){
#   print(bwdf$target[i])
#   print(bwdf$totalScore[i])
#   mcols(gr)[paste0(bwdf$target[i],"_signalNorm")]<-NULL#round(1e9*as.matrix(mcols(gr)[paste0(bwdf$target[i],"_signal")])/bwdf$totalScore[i],0)
# }

saveRDS(data.frame(gr),paste0(workDir,runName,"/custom/rds/",prefix,"_ChIPpeaks.results_annotated.RDS"))

res<-readRDS(paste0(workDir,runName,"/custom/rds/",prefix,"_ChIPpeaks.results_annotated.RDS"))


ss<-res %>% filter(chrRegionType %in% c("arm","tip"),
                   seqnames %in% seqnames(Celegans)[1:5],
                   group=="Double mutant")

peakCols<-grep("_peaks",colnames(ss))
binMat<-ss[,peakCols]>0
row.names(binMat)<-ss$gene_id
binMat
colors <- brewer.pal(4, "Accent")
named_colors <- setNames(colors,colnames(binMat))
fit<-euler(binMat)
p1<-plot(fit,quantities=T, fills = list(fill = named_colors, alpha = 0.7))
p1


ss<-res %>% filter(chrRegionType %in% c("arm","tip"),
                   seqnames %in% seqnames(Celegans)[1:5],
                   peakTypeNumber==3)
table(ss$group)

obsCounts<-data.frame(ss) %>% group_by(group) %>%
  summarize(count = n())
tt<-data.frame(ss) %>% t_test(log2FoldChange~group,var.equal=F) %>%
  adjust_pvalue(method="fdr") %>% p_format(p.adj, new.col=T,accuracy=1e-32)
ylimits<-c(-1.5,2)

p2<-ggplot(ss,aes(x=group,y=log2FoldChange)) +
  geom_boxplot(aes(fill=group),outlier.shape=NA) +
  coord_cartesian(ylim=ylimits) +
  geom_hline(yintercept=0,linetype="dashed",color="grey40")+
  geom_text(data=obsCounts,aes(label=count,y=ylimits[1]*0.95),color="blue",angle=0,size=3) +
  stat_pvalue_manual(data=tt,label="p.adj.format",y.position=ylimits[2]*0.77,step.increase=0.035,
                     angle=0, hide.ns=T,tip.length=0.01,label.size=3) +
  labs(title=paste0(length(unique(ss$gene_id)),
                    " genes on autosomal arms<br>with HPL-2/LIN-61/H3K9me2 ChIP peaks"),
       y="log<sub>2</sub>FC")+
  theme(legend.position = "none", axis.text.x = element_text(angle = 90, hjust = 1),
        axis.title.x=element_blank())
p2

p<-ggarrange(p1,p2,ncol=2)
p
ggsave(filename = paste0(workDir, runName, "/custom/finalFigures/integrationWithChIP/lfcByGroup_autosomalArmTripleChIPpeaks",
                         "_padj",padjVal,"_lfc",lfcVal,".pdf"),
       plot = p, width = 15, height = 12,units="cm")


# create quantiles
ss
for(i in 1:nrow(bwdf)){
  print(bwdf$target[i])
  print(bwdf$totalScore[i])
  chipValues<-as.matrix(ss[paste0(bwdf$target[i],"_signal")])
  ss[,paste0(bwdf$target[i],"_quantile")]<-factor(cut(chipValues,breaks=quantile(chipValues,probs=c(0,0.33,0.66,1)),labels=F,include.lowest=T))
}

ylimits<-c(-1.5,2.5)
obsCounts<-data.frame(ss) %>% group_by(group,hpl2_quantile) %>%
  summarize(count = n())
tt<-data.frame(ss) %>% group_by(hpl2_quantile) %>% wilcox_test(log2FoldChange~group) %>%
  adjust_pvalue(method="fdr") %>% p_format(p.adj, new.col=T,accuracy=1e-32) %>%
  mutate(y.position=rep(seq(ylimits[2]*0.5,ylimits[2]*0.95,ylimits[2]*0.45/5),3))


p3<-ggplot(ss,aes(x=group,y=log2FoldChange)) +
  geom_boxplot(aes(fill=group),outlier.shape=NA) +
  coord_cartesian(ylim=ylimits) +
  facet_wrap(.~hpl2_quantile)+
  geom_hline(yintercept=0,linetype="dashed",color="grey40")+
  geom_text(data=obsCounts,aes(label=count,y=ylimits[1]*0.95),color="blue",angle=0,size=3) +
  stat_pvalue_manual(data=tt,label="p.adj.format",
                     angle=0, hide.ns=T,tip.length=0.01,label.size=3) +
  labs(title=paste0(length(unique(ss$gene_id)),
                    " genes on autosomal arms with HPL-2/LIN-61/H3K9me2 ChIP peaks<br>grouped by HPL-2 ChIP signal quantile"),
       y="log<sub>2</sub>FC")+
  theme(legend.position = "none",
        axis.title.x=element_blank(),axis.text.x=element_blank())
p3


ylimits<-c(-1.5,2.5)
obsCounts<-data.frame(ss) %>% group_by(group,lin61_quantile) %>%
  summarize(count = n())
tt<-data.frame(ss) %>% group_by(lin61_quantile) %>%  wilcox_test(log2FoldChange~group) %>%
  adjust_pvalue(method="fdr") %>% p_format(p.adj, new.col=T,accuracy=1e-32) %>%
  mutate(y.position=rep(seq(ylimits[2]*0.5,ylimits[2]*0.95,ylimits[2]*0.45/5),3))

p4<-ggplot(ss,aes(x=group,y=log2FoldChange)) +
  geom_boxplot(aes(fill=group),outlier.shape=NA) +
  coord_cartesian(ylim=ylimits) +
  facet_wrap(.~lin61_quantile)+
  geom_hline(yintercept=0,linetype="dashed",color="grey40")+
  geom_text(data=obsCounts,aes(label=count,y=ylimits[1]*0.95),color="blue",angle=0,size=3) +
  stat_pvalue_manual(data=tt,label="p.adj.format",
                     angle=0, hide.ns=T,tip.length=0.01,label.size=3) +
  labs(title=paste0(length(unique(ss$gene_id)),
                    " genes on autosomal arms with HPL-2/LIN-61/H3K9me2 ChIP peaks<br>grouped by LIN-61 ChIP signal quantile"),
       y="log<sub>2</sub>FC")+
  theme(legend.position = "none",
        axis.title.x=element_blank(),axis.text.x=element_blank())
p4


ylimits<-c(-1.5,2.5)
obsCounts<-data.frame(ss) %>% group_by(group,H3K9me2_quantile) %>%
  summarize(count = n())
tt<-data.frame(ss) %>% group_by(H3K9me2_quantile) %>%  wilcox_test(log2FoldChange~group) %>%
  adjust_pvalue(method="fdr") %>% p_format(p.adj, new.col=T,accuracy=1e-32) %>%
  mutate(y.position=rep(seq(ylimits[2]*0.5,ylimits[2]*0.95,ylimits[2]*0.45/5),3))

p5<-ggplot(ss,aes(x=group,y=log2FoldChange)) +
  geom_boxplot(aes(fill=group),outlier.shape=NA) +
  coord_cartesian(ylim=ylimits) +
  facet_wrap(.~H3K9me2_quantile)+
  geom_hline(yintercept=0,linetype="dashed",color="grey40")+
  geom_text(data=obsCounts,aes(label=count,y=ylimits[1]*0.95),color="blue",angle=0,size=3) +
  stat_pvalue_manual(data=tt,label="p.adj.format",
                     angle=0, hide.ns=T,tip.length=0.01,label.size=3) +
  labs(title=paste0(length(unique(ss$gene_id)),
                    " genes on autosomal arms with HPL-2/LIN-61/H3K9me2 ChIP peaks<br>grouped by H3K9me2 ChIP signal quantile"),
       y="log<sub>2</sub>FC")+
  theme(legend.position = "none", axis.text.x = element_text(angle = 90, hjust = 1),
        axis.title.x=element_blank())
p5

p<-ggarrange(p3,p4,p5,nrow=3,heights=c(1.3,1.3,2))
p
ggsave(filename = paste0(workDir, runName, "/custom/finalFigures/integrationWithChIP/lfcByGroup_autosomalArmTripleChIPpeaks",
                         "_padj",padjVal,"_lfc",lfcVal,"_quantile_wilcoxTest.pdf"),
       plot = p, width = 15, height = 30,units="cm")

