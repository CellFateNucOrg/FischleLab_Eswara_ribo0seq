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
runName="/diff_abund_2_canonical_noRRnoSP"

contrasts<-read.csv(paste0(workDir,"/contrasts.csv"),sep=",",header=T)
prefix="ribo0_canonical_geneset_all"

setwd(workDir)

dir.create(paste0(workDir,runName,"/custom/plots"), showWarnings = FALSE, recursive = TRUE)


### functions ------
tableToGranges<-function(results,sort=TRUE){
  gr <- makeGRangesFromDataFrame(results, keep.extra.columns=TRUE,
                                 seqnames.field="seqnames",
                                 start.field="start",
                                 end.field="end",
                                 strand.field="strand")
  if(sort){
    gr<-sort(gr)
  }
  return(gr)
}



## interactive Volcano plots -------
results<-readRDS(paste0(workDir,runName,"/custom/rds/",prefix,".results_annotated.RDS"))
dir.create(paste0(workDir,runName,"/custom/plots/volcano"), showWarnings = FALSE, recursive = TRUE)
rockman<-readRDS(paste0(serverPath,"/publicData/Various/chrRegions_Rockman2009.RDS"))
gr<-tableToGranges(results,sort=FALSE)
seqlevelsStyle(gr)<-"UCSC"
ol<-findOverlaps(resize(gr,width=1,fix="start"),rockman,ignore.strand=T)
gr$chrRegionType[queryHits(ol)] <- paste0(rockman$chrRegionType[subjectHits(ol)])
gr$chrRegion[queryHits(ol)] <- paste0(seqnames(rockman)[subjectHits(ol)],"_",rockman$chrRegionType[subjectHits(ol)])
res<-as.data.frame(gr)
contrasts

lfcVal=1
padjVal=0.05

# free scale
for(i in 1:length(contrasts$id)){
  tmp<-res[res$shortId==contrasts$shortId[i],c("padj","log2FoldChange","gene_name","gene_biotype","chrRegion","chrRegionType")]
  tmp$mlog10padj<--log10(tmp$padj)
  tmp$significant<-"NS"
  tmp$significant[tmp$log2FoldChange>lfcVal & tmp$padj<padjVal]<- "Up"
  tmp$significant[tmp$log2FoldChange< -lfcVal & tmp$padj<padjVal]<- "Down"

  tmp$significant<-factor(tmp$significant,levels=c("Up","Down","NS"))


  volcano_plot <- ggplot(tmp, aes(x = log2FoldChange, y = -log10(padj),
                                  color=significant, shape = chrRegionType,
                                  text = paste("Gene:", gene_name,
                                               "<br>Biotype:", gene_biotype,
                                               "<br>chrRegion:", chrRegion,
                                               "<br>log2FC:", round(log2FoldChange, 2),
                                               "<br>padj:", signif(padj, 3)))) +
    geom_point(size=1.5,alpha=0.4) +
    scale_color_manual(values = c("Up" = "#E41A1C",
                                  "Down" = "#377EB8",
                                  "NS" = "lightgrey"),
                       guide = "none") +
    theme_minimal() + ggtitle(contrasts$id[i]) + xlab("Log2 Fold Change")

  volcano_plot1<-volcano_plot + geom_text_repel(
    data = subset(tmp, padj < 0.05 & abs(log2FoldChange) > 0.5),
    aes(label = gene_name),
    box.padding = 0,
    #max.overlaps = Inf,
    size = 2,
    force = 0.5, force_pull=1
  )
  volcano_plot1
  ggsave(filename = paste0(workDir,runName,"/custom/plots/volcano/volcano_",contrasts$id[i],"_free.png"),
         plot = volcano_plot1, width = 8, height = 6, dpi = 300,bg="white")
  volcano_plot
  interactive_volcano <- ggplotly(volcano_plot,tooltip="text")
  saveWidget(interactive_volcano, file = paste0(workDir,runName,"/custom/plots/volcano/volcano_",contrasts$id[i],"_free.html"))
}

# zoom LFC -5 to 5
for(i in 1:length(contrasts$id)){
  tmp<-res[res$shortId==contrasts$shortId[i],c("padj","log2FoldChange","gene_name","gene_biotype","chrRegion","chrRegionType")]
  tmp$mlog10padj<--log10(tmp$padj)
  tmp$significant<-"NS"
  tmp$significant[tmp$log2FoldChange>lfcVal & tmp$padj<padjVal]<- "Up"
  tmp$significant[tmp$log2FoldChange< -lfcVal & tmp$padj<padjVal]<- "Down"

  tmp$significant<-factor(tmp$significant,levels=c("Up","Down","NS"))


  volcano_plot <- ggplot(tmp, aes(x = log2FoldChange, y = -log10(padj),
                                  color=significant, shape = chrRegionType,
                                  text = paste("Gene:", gene_name,
                                               "<br>Biotype:", gene_biotype,
                                               "<br>chrRegion:", chrRegion,
                                               "<br>log2FC:", round(log2FoldChange, 2),
                                               "<br>padj:", signif(padj, 3)))) +
    geom_point(size=1.5,alpha=0.4) +
    scale_color_manual(values = c("Up" = "#E41A1C",
                                  "Down" = "#377EB8",
                                  "NS" = "lightgrey"),
                       guide = "none") +
    theme_minimal() + ggtitle(contrasts$id[i]) + xlab("Log2 Fold Change") +
    coord_cartesian(xlim=c(-5,5))

  volcano_plot1<-volcano_plot + geom_text_repel(
    data = subset(tmp, padj < 0.05 & abs(log2FoldChange) > 0.5),
    aes(label = gene_name),
    box.padding = 0,
    #max.overlaps = Inf,
    size = 2,
    force = 0.5, force_pull=1
  )
  volcano_plot1
  ggsave(filename = paste0(workDir,runName,"/custom/plots/volcano/volcano_",contrasts$id[i],"_-5to5.png"),
         plot = volcano_plot1, width = 8, height = 6, dpi = 300,bg="white")
  volcano_plot
  interactive_volcano <- ggplotly(volcano_plot,tooltip="text")
  saveWidget(interactive_volcano, file = paste0(workDir, runName,"/custom/plots/volcano/volcano_",contrasts$id[i],"_-5to5.html"))
}


## Up down by chromosome  -----
lfcVal=1
padjVal=0.05


dir.create(paste0(workDir, runName, "/custom/plots/byChromosome"), showWarnings = FALSE, recursive = TRUE)
nuclear<-seqlevels(Celegans)[1:6]
resA<-res[res$seqnames %in% nuclear,]
resA$padj[is.na(resA$padj)]<-1

resA<-resA[resA$padj<padjVal & abs(resA$log2FoldChange)>lfcVal,]
obsCounts<-data.frame(resA) %>% group_by(group,seqnames) %>%
  summarize(count = n())
wilcoxt<-data.frame(resA) %>% group_by(group) %>% wilcox_test(log2FoldChange~seqnames,ref.group="all") %>%
  adjust_pvalue(method="fdr") %>% p_format(new.col=T,accuracy=1e-32)
ylimits<-c(-3,3)
p<-ggplot(resA,aes(x=seqnames,y=log2FoldChange)) +
  geom_boxplot(aes(fill=seqnames),outlier.shape=NA) +
  facet_wrap(~group, scales="free") +
  coord_cartesian(ylim=ylimits) +
  geom_hline(yintercept=0,linetype="dashed",color="grey40") +
  scale_fill_brewer(palette="Blues")+
  geom_text(data=obsCounts,aes(label=count,y=ylimits[1]*0.95),color="blue",angle=0,size=3) +
  stat_pvalue_manual(data=wilcoxt,label="p.adj.signif",y.position=ylimits[2]*0.95,
                     remove.bracket = T,angle=0,hide.ns=T) +
  ggtitle(paste0("Log2FoldChange for significant genes (padj<",padjVal," |LFC|>",lfcVal,")"))+
  theme(legend.position = "none")
p
ggsave(filename = paste0(workDir, runName, "/custom/plots/byChromosome/lfcByChr_",
                         contrasts$id[i], "_padj",padjVal,"_lfc",lfcVal,".png"),
       plot = p, width = 11, height = 8, dpi = 300,bg="white")

resA$upVdown<-NA
resA$upVdown[resA$log2FoldChange>lfcVal & resA$padj<padjVal]<-"up"
resA$upVdown[resA$log2FoldChange<lfcVal & resA$padj<padjVal]<-"down"
resA$upVdown<-factor(resA$upVdown, levels=c("up","down"))

obsCounts<-data.frame(resA) %>% group_by(group,seqnames,upVdown) %>%
  summarize(count = n())
obsCounts$y<-ifelse(obsCounts$upVdown=="up",0.95,0.05)
p<-ggplot(resA,aes(x=seqnames,fill=upVdown)) +
  geom_bar(position="fill") +
  facet_wrap(~group) +
  geom_hline(yintercept=0.5,linetype="dashed",color="black") +
  geom_text(data=obsCounts,aes(label=count,y=y),color="black",angle=0, size=3) +
  ggtitle(paste0("Fraction of up/down significant genes (padj<",padjVal," |LFC|>",lfcVal,")"))+
  scale_fill_brewer(palette = "Accent")
p
ggsave(filename = paste0(workDir, runName, "/custom/plots/byChromosome/countsByChr_",
                         contrasts$id[i],"_padj",padjVal,"_lfc",lfcVal,".png"),
       plot = p, width = 11, height = 8, dpi = 300,bg="white")




## Up down chromosome regions -----
dir.create(paste0(workDir, runName,"/custom/plots/byChrRegion"), showWarnings = FALSE, recursive = TRUE)
autosomes<-seqlevels(Celegans)[1:5]
resA<-res[res$seqnames %in% autosomes,]
resA$padj[is.na(resA$padj)]<-1
resA$chrRegionType<-factor(resA$chrRegionType, levels=c("tip","arm","center"))

resA<-resA[resA$padj<padjVal & abs(resA$log2FoldChange)>lfcVal,]
obsCounts<-data.frame(resA) %>% group_by(group,chrRegionType) %>%
  summarize(count = n())
wilcoxt<-data.frame(resA) %>% group_by(group) %>% wilcox_test(log2FoldChange~chrRegionType,ref.group="all") %>%
  adjust_pvalue(method="fdr") %>% p_format(new.col=T,accuracy=1e-32)
ylimits<-c(-3,3)
p<-ggplot(resA,aes(x=chrRegionType,y=log2FoldChange)) +
  geom_boxplot(aes(fill=chrRegionType),outlier.shape=NA) +
  facet_wrap(~group, scales="free") +
  coord_cartesian(ylim=ylimits) +
  geom_hline(yintercept=0,linetype="dashed",color="grey") +
  scale_fill_manual(values=c("tip"="blue","arm"="royalblue1","center"="lightblue"))+
  geom_text(data=obsCounts,aes(label=count,y=ylimits[1]*0.95),color="blue",angle=0) +
  stat_pvalue_manual(data=wilcoxt,label="p.adj.signif",y.position=ylimits[2]*0.95,
                     remove.bracket = T,angle=0,hide.ns=T) +
  ggtitle(paste0("Log2FoldChange for significant autosomal genes (padj<",padjVal," |LFC|>",lfcVal,")"))
ggsave(filename = paste0(workDir, runName, "/custom/plots/byChrRegion/lfcByChrRegion_",
                         contrasts$id[i], "_padj",padjVal,"_lfc",lfcVal,".png"),
       plot = p, width = 11, height = 8, dpi = 300,bg="white")

resA$upVdown<-NA
resA$upVdown[resA$log2FoldChange>lfcVal & resA$padj<padjVal]<-"up"
resA$upVdown[resA$log2FoldChange<lfcVal & resA$padj<padjVal]<-"down"
resA$upVdown<-factor(resA$upVdown, levels=c("up","down"))

obsCounts<-data.frame(resA) %>% group_by(group,chrRegionType,upVdown) %>%
  summarize(count = n())
obsCounts$y<-ifelse(obsCounts$upVdown=="up",0.95,0.05)
p<-ggplot(resA,aes(x=chrRegionType,fill=upVdown)) +
  geom_bar(position="fill") +
  facet_wrap(~group) +
  geom_hline(yintercept=0.5,linetype="dashed",color="black") +
  geom_text(data=obsCounts,aes(label=count,y=y),color="black",angle=0) +
  ggtitle(paste0("Fraction of up/down significant autosomal genes (padj<",padjVal," |LFC|>",lfcVal,")"))+
  scale_fill_brewer(palette = "Accent")

ggsave(filename = paste0(workDir, runName, "/custom/plots/byChrRegion/countsByChrRegion_",
                         contrasts$id[i],"_padj",padjVal,"_lfc",lfcVal,".png"),
       plot = p, width = 10, height = 8, dpi = 300,bg="white")





## Heatmap of significant genes -----
results<-readRDS(paste0(workDir,runName,"/custom/rds/",prefix,".results_annotated.RDS"))
dir.create(paste0(workDir, runName, "/custom/plots/heatmaps"), showWarnings = FALSE, recursive = TRUE)
n2contrasts<-grep("vs_N2$",levels(results$group),value=T)
lin61contrasts<-grep("vs_HPL2GFP__lin61$",levels(results$group),value=T)

res_n2<-results[results$group %in% n2contrasts,]
res_lin61<-results[results$group %in% lin61contrasts,]

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


getLFCMat<-function(results, numSamplesSignificant=1, lfcVal=0, padjVal=0.05){
  mat_padj<-gatherResults(results,valueColumn="padj")
  mat_padj[is.na(mat_padj)]<-1
  mat_lfc<-gatherResults(results,valueColumn="log2FoldChange")

  sig<-(mat_padj[,2:ncol(mat_padj)]<padjVal & abs(mat_lfc[,2:ncol(mat_lfc)])>lfcVal)
  sigGenes<-mat_padj$gene_id[rowSums(sig)>numSamplesSignificant]

  mat_lfc<-mat_lfc[mat_lfc$gene_id %in% sigGenes,]
  rownames(mat_lfc)<-mat_lfc$gene_id
  mat_lfc<-mat_lfc[,2:ncol(mat_lfc)]
  mat_lfc<-as.matrix(mat_lfc)
  colnames(mat_lfc)<-gsub("log2FoldChange_","",colnames(mat_lfc))
  return(mat_lfc)
}




numSamplesSignificant=3
lfcVal=1
padjVal=0.05
mat_lfc<-getLFCmat(res_n2,numSamplesSignificant, lfcVal, padjVal)

png(paste0(workDir, runName,"/custom/plots/heatmaps/hclust_heatmap_sigSamples",
           numSamplesSignificant,"_padj",padjVal,"_lfc",lfcVal,"_n2Contrasts.png"),
    width=19,height=29,units="cm",res=150)
ht<-Heatmap(mat_lfc,show_row_names=F, cluster_columns=T, cluster_rows=T, show_row_dend = F,
            column_names_rot=75,column_title=paste0("log2FC of ",nrow(mat_lfc),
                                                    " genes singificant in >=",numSamplesSignificant,
                                                    " samples (padj<",padjVal," |LFC|>",lfcVal,")"),
            col=circlize::colorRamp2(breaks = c(-4, 0, 4), colors = c("blue", "white", "red")))
ht<-draw(ht)
dev.off()


mat_lfc<-getLFCmat(res_lin61, numSamplesSignificant, lfcVal, padjVal)

png(paste0(workDir, runName, "/custom/plots/heatmaps/hclust_heatmap_sigSamples",
           numSamplesSignificant,"_padj",padjVal,"_lfc",lfcVal,"_lin61Contrasts.png"),
    width=19,height=29,units="cm",res=150)
ht<-Heatmap(mat_lfc,show_row_names=F, cluster_columns=T, cluster_rows=T, show_row_dend = F,
            column_names_rot=75, column_title=paste0("log2FC of ",nrow(mat_lfc),
                                                     " genes singificant in >=",numSamplesSignificant,
                                                     " samples (padj<",padjVal," |LFC|>",lfcVal,")"),
            col=circlize::colorRamp2(breaks = c(-4, 0, 4), colors = c("blue", "white", "red")))
ht<-draw(ht)
dev.off()




mat_lfc<-getLFCmat(results, numSamplesSignificant, lfcVal, padjVal)

png(paste0(workDir,runName,"/custom/plots/heatmaps/hclust_heatmap_sigSamples",
           numSamplesSignificant,"_padj",padjVal,"_lfc",lfcVal,"_allContrasts.png"),
    width=19,height=29,units="cm",res=150)
ht<-Heatmap(mat_lfc,show_row_names=F, cluster_columns=T, cluster_rows=T, show_row_dend = F,
            column_names_rot=90, column_title=paste0("log2FC of ",nrow(mat_lfc),
                                                     " genes singificant in >=",numSamplesSignificant,
                                                     " samples (padj<",padjVal," |LFC|>",lfcVal,")"),
            col=circlize::colorRamp2(breaks = c(-4, 0, 4), colors = c("blue", "white", "red")))
ht<-draw(ht)
dev.off()

