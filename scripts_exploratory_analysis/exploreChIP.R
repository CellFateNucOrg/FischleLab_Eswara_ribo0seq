library(GenomicRanges)
library(rtracklayer)
library(ggplot2)
library(BSgenome.Celegans.UCSC.ce11)
library(ggpubr)
library(rstatix)
library(dplyr)
library(RColorBrewer)
library(tidyr)

theme_set(
  theme_classic()+
    theme(panel.grid.major = element_blank(),
          panel.grid.minor = element_blank(),
          axis.title.y=ggtext::element_markdown(size=9),
          axis.title.x=ggtext::element_markdown(size=9),
          title=ggtext::element_markdown(size=9)
    )
)
options(tibble.width=Inf)

serverPath="/Volumes/external.data/MeisterLab"
#serverPath="Z:/MeisterLab"

workDir=paste0(serverPath,"/FischleLab_KarthikEswara/ribo0seq")
runName="/diff_abund_3_canonical_noRRnoSP_moreSeq"

source(paste0(workDir,"/scripts_exploratory_analysis/functions_exploratory_analysis.R"))
contrasts<-read.csv(paste0(workDir,"/contrasts_subset.csv"),sep=",",header=T)

prefix="ribo0_canonical_geneset_all"

setwd(workDir)

dir.create(paste0(workDir,runName,"/custom/exploreChIP"),showWarnings = F,recursive = T)

lfcVal=0.5
padjVal=0.05

## deeptools ------
## bed files for deeptools
results<-readRDS(paste0(workDir,runName,"/custom/rds/",prefix,"_allArmGenes.results_annotated.RDS"))
dir.create(paste0(workDir,runName,"/custom/finalFigures/exploreChIP/bed"),showWarnings = F,recursive = T)
for(g in levels(results$group)){
 tmp<-results[results$group==g,]
 gr<-GRanges(tmp)
 gr$name<-tmp$gene_id
 gr$score<-tmp$log2FoldChange
 gr$padj[is.na(gr$padj)]<-1
 forBed<-gr[gr$padj<padjVal & gr$log2FoldChange>lfcVal]
 export.bed(forBed[order(forBed$log2FoldChange,decreasing =T)],paste0(workDir,"/",runName,"/custom/exploreChIP/bed/",gsub(" ","_",g),"__autosomalArmsUp_",prefix,".bed"))
 forBed<-gr[gr$padj<padjVal & gr$log2FoldChange < (-lfcVal)]
 export.bed(forBed[order(forBed$log2FoldChange,decreasing = F)],paste0(workDir,"/",runName,"/custom/exploreChIP/bed/",gsub(" ","_",g),"__autosomalArmsDown_",prefix,".bed"))
}

bwList<-read.csv(paste0(workDir,"/publicChIPseq.txt"),header=F)
if(ncol(bwList)==1){
  target<-sapply(strsplit(sapply(bwList,basename),"_"),"[[",1)
  id<-sapply(strsplit(sapply(bwList,basename),"_"),"[[",2)
  df<-data.frame(target=target,id=id,bigwig=unlist(bwList),stringsAsFactors = F)
  bworder<-c(7,9,5,6,8,2,1,3,4)
  write.table(df[bworder,],paste0(workDir,"/",runName,"/custom/exploreChIP/publicChIPseq.tsv"),
              row.names=F,col.names=F,quote=F,sep="\t")
}


## boxplots - peak overlap with 10kb region -----
### all genes
results<-readRDS(paste0(workDir,runName,"/custom/rds/",prefix,".results_annotated.RDS"))
results$longId<-droplevels(results$group)
results$shortId<-droplevels(results$shortId)
results$group<-factor(results$shortId,levels=contrasts$shortId,labels=contrasts$prettyName)
results<-results[results$group %in% contrasts$prettyName,]
results<-results[results$seqnames!="MtDNA",]

gr<-tableToGranges(results)
triplePeaks<-readRDS(paste0(workDir,"/../infoFromKarthik/ce11_triple_peaks_db.rds"))
ol<-findOverlaps(resize(gr,width=1,fix="center"),triplePeaks,ignore.strand=T)
gr$triplePeakClass<-NA
gr$triplePeakClass[queryHits(ol)]<-triplePeaks$triplePeaks_WF_density_class[subjectHits(ol)]


seqlevelsStyle(gr)<-"UCSC"
rockman<-readRDS(paste0(serverPath,"/publicData/Various/chrRegions_Rockman2009.RDS"))
ol<-findOverlaps(resize(gr,width=1,fix="start"),rockman,ignore.strand=T)
gr$chrRegionType[queryHits(ol)] <- paste0(rockman$chrRegionType[subjectHits(ol)])
gr$chrRegion[queryHits(ol)] <- paste0(seqnames(rockman)[subjectHits(ol)],"_",rockman$chrRegionType[subjectHits(ol)])



obsCounts<-data.frame(gr) %>% group_by(chrRegionType,group,triplePeakClass) %>%
  summarize(count = n())
tt<-data.frame(gr) %>% group_by(chrRegionType, group) %>%
  t_test(log2FoldChange~triplePeakClass, ref.group = "all", var.equal=F) %>%
  adjust_pvalue(method="fdr")  %>% rstatix::p_format(p.adj,new.col=T,accuracy=1e-32)
ylimits<-c(-1.5,1.5)

p<-ggplot(as.data.frame(gr),aes(x=triplePeakClass,y=log2FoldChange,fill=triplePeakClass)) +
         geom_boxplot(outlier.shape=NA)+
  facet_grid(chrRegionType~group) +
  coord_cartesian(ylim=ylimits) +
  geom_hline(yintercept=0,linetype="dashed",color="grey")+
  geom_text(data=obsCounts,aes(label=count,y=ylimits[1]*0.95),color="blue",angle=0,size=3) +
  stat_pvalue_manual(data=tt,label="p.adj.format",y.position=ylimits[2]*0.95,
                     remove.bracket=T,angle=0, hide.ns="p.adj",  label.size=3) +
  labs(title=paste0(length(unique(gr$gene_id)),
                    " (all genes) overlapping 10kb triple peaks class"),
       y="log<sub>2</sub>FC")+
  theme(legend.position = "none", axis.text.x = element_text(angle = 90, hjust = 1),
        axis.title.x=element_blank())
p
ggsave(paste0(workDir,runName,"/custom/exploreChIP/allGenes_triplePeaks10kb_boxplot.pdf"),
       p,width=8,height=8)


### gene subset
res<-readRDS(paste0(workDir,runName,"/custom/rds/",prefix,"_subset.results_annotated.RDS"))
table(res$group)
gr<-tableToGranges(res)
triplePeaks<-readRDS(paste0(workDir,"/../infoFromKarthik/ce11_triple_peaks_db.rds"))
seqinfo(triplePeaks)
seqlevels(triplePeaks)<-paste0("chr",seqlevels(triplePeaks))
seqlevels(triplePeaks)<-c(seqlevels(triplePeaks),"chrM")


ol<-findOverlaps(resize(gr,width=1,fix="center"),triplePeaks,ignore.strand=T)
gr$triplePeakClass<-NA
gr$triplePeakClass[queryHits(ol)]<-triplePeaks$triplePeaks_WF_density_class[subjectHits(ol)]


obsCounts<-data.frame(gr) %>% group_by(group,triplePeakClass) %>%
  summarize(count = n())
tt<-data.frame(gr) %>% group_by(group) %>%
  t_test(log2FoldChange~triplePeakClass, ref.group = "all", var.equal=F) %>%
  adjust_pvalue(method="fdr")  %>% rstatix::p_format(p.adj,new.col=T,accuracy=1e-32)
ylimits<-c(0,2.5)

p<-ggplot(as.data.frame(gr),aes(x=triplePeakClass,y=log2FoldChange,fill=triplePeakClass)) +
  geom_boxplot(outlier.shape=NA)+
  facet_grid(.~group) +
  coord_cartesian(ylim=ylimits) +
  geom_hline(yintercept=0,linetype="dashed",color="grey")+
  geom_text(data=obsCounts,aes(label=count,y=ylimits[1]*0.95+0.1),color="blue",angle=0,size=3) +
  stat_pvalue_manual(data=tt,label="p.adj.format",y.position=ylimits[2]*0.95,
                     remove.bracket=T,angle=0, hide.ns="p.adj",  label.size=3) +
  labs(title=paste0(length(unique(gr$gene_id)),
                    " autosomal arm significant genes (padj<",padjVal," LFC>",lfcVal,")<br>overlapping 10kb triple peaks class"),
       y="log<sub>2</sub>FC")+
  theme(legend.position = "none", axis.text.x = element_text(angle = 90, hjust = 1),
        axis.title.x=element_blank())
p
ggsave(paste0(workDir,runName,"/custom/exploreChIP/autosomalArms_padj",padjVal,"_lfc",
              lfcVal,"_triplePeaks10kb_boxplot.pdf"),
       p,width=8,height=4)

## count Peak overlap ----
bedFiles<-list.files(paste0(workDir,"/../GSE271919_ChIPseq/bed/"),pattern=".*GSE.*\\.bed")
targets<-sapply(strsplit(bedFiles,"_"),"[[",1)
chosen<-c(7,9,5)
beddf<-data.frame(bedFiles=paste0(workDir,"/../GSE271919_ChIPseq/bed/",bedFiles[chosen]),targets=targets[chosen])

### all genes
results<-readRDS(paste0(workDir,runName,"/custom/rds/",prefix,".results_annotated.RDS"))
results$padj[is.na(results$padj)]<-1
results$longId<-droplevels(results$group)
results$shortId<-droplevels(results$shortId)
results$group<-factor(results$shortId,levels=contrasts$shortId,labels=contrasts$prettyName)
results<-results[results$group %in% contrasts$prettyName,]
results<-results[results$seqnames!="MtDNA",]

results<-results[results$group=="Double mutant",]

gr<-tableToGranges(results)
head(gr)
#head(resize(gr,width=width(gr)+100,fix="end"))

i=1
for(i in 1:nrow(beddf)){
  bed<-import(beddf$bedFiles[i])
  cname<-paste0(beddf$target[i],"_peaks")
  mcols(gr)[cname]<-countOverlaps(gr,bed)
}


seqlevelsStyle(gr)<-"UCSC"
rockman<-readRDS(paste0(serverPath,"/publicData/Various/chrRegions_Rockman2009.RDS"))
ol<-findOverlaps(resize(gr,width=1,fix="start"),rockman,ignore.strand=T)
gr$chrRegionType[queryHits(ol)] <- paste0(rockman$chrRegionType[subjectHits(ol)])
gr$chrRegion[queryHits(ol)] <- paste0(seqnames(rockman)[subjectHits(ol)],"_",rockman$chrRegionType[subjectHits(ol)])

df<-tidyr::pivot_longer(data.frame(gr), cols=c("hpl2_peaks","H3K9me2_peaks","lin61_peaks"),
             names_to="peakSet",values_to="peakCount")

head(df)

p1<-ggplot(df,aes(x=peakCount,fill=peakSet)) +
  geom_histogram()+
  facet_grid(peakSet~chrRegionType) +
  labs(title=paste0(length(unique(df$gene_id)),
                    " (all genes) overlapping bed file peaks per gene"),
       y="log<sub>2</sub>FC")+
  theme(legend.position = "none", axis.text.x = element_text(angle = 90, hjust = 1),
        axis.title.x=element_blank())
p1

p2<-ggplot(df[df$peakCount>0,],aes(x=peakCount,fill=peakSet)) +
  geom_histogram()+
  facet_grid(peakSet~chrRegionType) +
  labs(title=paste0(length(unique(df$gene_id)),
                    " (all genes) overlapping bed file peaks per gene (>0)"),
       y="log<sub>2</sub>FC")+
  theme(legend.position = "none", axis.text.x = element_text(angle = 90, hjust = 1),
        axis.title.x=element_blank())
p2

p<-ggarrange(p1,p2,nrow=2)
ggsave(paste0(workDir,runName,"/custom/exploreChIP/allGenes_bedPeakOverlap_histogram.pdf"),
       p,width=8,height=12)

#### overlap of peak types ----
head(gr)
peakCols<-grep("_peaks$",colnames(mcols(gr)))
doublePeakCols<-peakCols[c(1,3)]

gr$peakTypeNum<-rowSums(data.frame(mcols(gr)[,peakCols])>0)
gr$peakTypeNum_double<-rowSums(data.frame(mcols(gr)[,doublePeakCols])>0)

df<-data.frame(gr)
df$peakTypeNum<-factor(df$peakTypeNum,levels=c(3:0))
df$peakTypeNum_double<-factor(df$peakTypeNum_double,levels=c(2:0))
p1<-ggplot(df) +
  geom_bar(aes(x=chrRegionType,group=peakTypeNum,
               fill=peakTypeNum),color="black")+
  labs(title=paste0("The number of peak types overlapping ",length(unique(df$gene_id)),
                    "<br>(all expressed) genes"))+
  theme(axis.text.x = element_text(angle = 90, hjust = 1),
        axis.title.x=element_blank()) +
  scale_fill_manual(values=c("black","grey40","grey80","white"))
p1

p1a<-ggplot(df) +
  geom_bar(aes(x=chrRegionType,group=peakTypeNum_double,
               fill=peakTypeNum_double),color="black")+
  labs(title=paste0("The number of hpl2/H3K9me2 peak types overlapping ",length(unique(df$gene_id)),
                    "<br>(all expressed) genes"))+
  theme(axis.text.x = element_text(angle = 90, hjust = 1),
        axis.title.x=element_blank()) +
  scale_fill_manual(values=c("black","grey60","white"))
p1a

### up genes
results<-readRDS(paste0(workDir,runName,"/custom/rds/",prefix,".results_annotated.RDS"))
results$padj[is.na(results$padj)]<-1
results<-results[results$shortId %in% contrasts$shortId,]
results$longId<-droplevels(results$group)
results$shortId<-droplevels(results$shortId)
results$group<-factor(results$shortId,levels=contrasts$shortId,labels=contrasts$prettyName)
results<-results[!(results$seqnames %in% c("MtDNA","X")),]
gr<-tableToGranges(results)
seqlevelsStyle(gr)<-"UCSC"
rockman<-readRDS(paste0(serverPath,"/publicData/Various/chrRegions_Rockman2009.RDS"))
ol<-findOverlaps(resize(gr,width=1,fix="start"),rockman,ignore.strand=T)
gr$chrRegionType[queryHits(ol)] <- paste0(rockman$chrRegionType[subjectHits(ol)])

gr<-gr[gr$chrRegionType %in% c("arm","tip")]

upGenes<-gr[gr$padj<padjVal & gr$log2FoldChange>lfcVal]
table(upGenes$group)
downGenes<-gr[gr$padj<padjVal & gr$log2FoldChange < (-lfcVal)]
table(downGenes$group)


for(i in 1:nrow(beddf)){
  bed<-import(beddf$bedFiles[i])
  seqlevelsStyle(bed)<-"UCSC"
  cname<-paste0(beddf$target[i],"_peaks")
  mcols(upGenes)[cname]<-countOverlaps(upGenes,bed)
  mcols(downGenes)[cname]<-countOverlaps(downGenes,bed)
}

table(upGenes$group)

peakCols<-grep("_peaks$",colnames(mcols(upGenes)))
doublePeakCols<-peakCols[c(1,3)]

upGenes$peakTypeNum<-rowSums(data.frame(mcols(upGenes)[,peakCols])>0)
downGenes$peakTypeNum<-rowSums(data.frame(mcols(downGenes)[,peakCols])>0)
upGenes$peakTypeNum_double<-rowSums(data.frame(mcols(upGenes)[,doublePeakCols])>0)
downGenes$peakTypeNum_double<-rowSums(data.frame(mcols(downGenes)[,doublePeakCols])>0)

df<-data.frame(upGenes)
df$peakTypeNum<-factor(df$peakTypeNum,levels=c(3:0))
df$peakTypeNum_double<-factor(df$peakTypeNum_double,levels=c(2:0))
p2<-ggplot(df) +
  geom_bar(aes(x=group,group=peakTypeNum,
               fill=peakTypeNum),color="black")+
  labs(title=paste0("The number of peak types overlapping ",length(unique(df$gene_id)),
                    " significantly<br>upregulated (padj<",padjVal," LFC>",lfcVal,
                    ") autosomal arm genes"))+
  theme(axis.text.x = element_text(angle = 90, hjust = 1),
        axis.title.x=element_blank()) +
  scale_fill_manual(values=c("black","grey40","grey80","white"))
p2

p2a<-ggplot(df) +
  geom_bar(aes(x=group,group=peakTypeNum_double,
               fill=peakTypeNum_double),color="black")+
  labs(title=paste0("The number of hpl2/H3K9me2 peak types overlapping ",length(unique(df$gene_id)),
                    " significantly<br>upregulated (padj<",padjVal," LFC>",lfcVal,
                    ") autosomal arm genes"))+
  theme(axis.text.x = element_text(angle = 90, hjust = 1),
        axis.title.x=element_blank()) +
  scale_fill_manual(values=c("black","grey60","white"))
p2a


df<-data.frame(downGenes)
df$peakTypeNum<-factor(df$peakTypeNum,levels=c(3:0))
df$peakTypeNum_double<-factor(df$peakTypeNum_double,levels=c(2:0))
p3<-ggplot(df) +
  geom_bar(aes(x=group,group=peakTypeNum,
               fill=peakTypeNum),color="black")+
  labs(title=paste0("The number of peak types overlapping ",length(unique(df$gene_id)),
                           " significantly<br>downregulated (padj<",padjVal," LFC< -",
                    lfcVal,") autosomal arm genes"))+
  theme(axis.text.x = element_text(angle = 90, hjust = 1),
        axis.title.x=element_blank()) +
  scale_fill_manual(values=c("black","grey40","grey80","white"))
p3

p3a<-ggplot(df) +
  geom_bar(aes(x=group,group=peakTypeNum_double,
               fill=peakTypeNum_double),color="black")+
  labs(title=paste0("The number of hpl2/H3K9me2 peak types overlapping ",length(unique(df$gene_id)),
                    " significantly<br>downregulated (padj<",padjVal," LFC< -",
                    lfcVal,") autosomal arm genes"))+
  theme(axis.text.x = element_text(angle = 90, hjust = 1),
        axis.title.x=element_blank()) +
  scale_fill_manual(values=c("black","grey60","white"))
p3a


p<-ggarrange(p1,p1a,p2,p2a,p3,p3a,nrow=3,ncol=2,widths=c(1,1),heights=c(0.9,1.1,1.1))
p
ggsave(paste0(workDir,runName,"/custom/exploreChIP/bedPeakTypeOverlap_barplot.pdf"),
       p,width=8,height=12)

