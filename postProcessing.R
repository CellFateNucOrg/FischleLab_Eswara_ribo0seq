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


### functions ------
annotate_results <- function(gtf,pattern="\\.deseq2\\.results\\.tsv",prefix=""){
  fileList<-list.files("./tables/differential",pattern=pattern, full.names=T)
  for(f in fileList){
    print(f)
    df <- read.delim(f, header=T)
    if(!("gene_id" %in% colnames(df))){
      df$gene_id<-row.names(df)
      df$gene_name<-row.names(df)
    }
    file_name <- gsub("\\.deseq2\\.results\\.tsv", "", basename(f))
    group_name <- gsub(paste0(prefix,"_"),"",file_name)
    print(group_name)
    df <- inner_join(df, data.frame(gtf), by=c("gene_id" = "gene_id"))
    df$group <- group_name
    write.table(df, paste0("./custom/txt/",file_name, ".deseq2.results_annotated.tsv"),
                row.names=F, quote=F, sep="\t", col.names=T)
  }
}


combine_results <- function(pattern="\\.deseq2\\.results_annotated\\.tsv"){
  fileList <- list.files("./custom/txt/", pattern=pattern, full.names=T)
  results <- do.call(rbind, lapply(fileList, read.delim, header=T))
  return(results)
}


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

#-----------------
contrasts<-read.csv("./contrasts.csv",sep=",",header=T)
prefix="ribo0_canonical_geneset_all"
setwd("/Volumes/external.data/MeisterLab/FischleLab_KarthikEswara/ribo0seq")
dir.create("custom/plots", showWarnings = FALSE, recursive = TRUE)
dir.create("custom/rds", showWarnings = FALSE, recursive = TRUE)
dir.create("custom/txt", showWarnings = FALSE, recursive = TRUE)
genomeVer<-"WS295"
gtf<-import(paste0("/Volumes/MeisterLab/publicData/genomes/",genomeVer,"/c_elegans.PRJNA13758.",genomeVer,".canonical_geneset.gtf"))
#gtf<-import("Z:/MeisterLab/publicData/genomes/WS295/c_elegans.PRJNA13758.WS295.canonical_geneset.gtf")
gtf <- gtf[gtf$type == "gene"]
mcols(gtf)<-mcols(gtf)[c("source","type","gene_id","gene_biotype","gene_name")]
gtf$source<-genomeVer
gtf<-sort(gtf)

annotate_results(gtf,pattern="\\.deseq2\\.results\\.tsv")

results<-combine_results(pattern="\\.deseq2\\.results_annotated\\.tsv")
results$seqnames<- factor(results$seqnames, levels=c("I", "II", "III", "IV", "V", "X", "MtDNA"))
results$shortId<-results$group
contrasts$shortId<- gsub("_","",contrasts$id)
#idx<-match(results$group,contrasts$shortId)
#results$contrast<-contrasts$id[idx]
results$shortId<-factor(results$group, levels=contrasts$shortId)
results$group<-factor(results$group, levels=contrasts$shortId, labels=contrasts$id)

if (file.exists(paste0("./custom/rds/",prefix,".results_annotated.RDS"))) {
  response <- readline(prompt = paste(prefix, "file already exists. Overwrite? [y/n]: "))
  if (tolower(response) == "y") {
    saveRDS(results, paste0("custom/rds/",prefix,".results_annotated.RDS"))
  } else {
    print("file not overwritten")
  }
} else {
  saveRDS(results, paste0("custom/rds/",prefix,".results_annotated.RDS"))
}



## Get counts of up/down by different thresholds-----
#contrasts<-read.csv("./contrasts.csv",sep=",",header=T)
# make table of up vs down
LFCthresholds<-c(0,0.5,1,1.5,2)
workDir<-"/Volumes/external.data/MeisterLab/FischleLab_KarthikEswara/ribo0seq/"
#tables<-list.files(path=paste0(workDir,"tables/differential"),pattern=".*deseq2\\.results\\.tsv")
contrastNames<-contrasts$id
tbl<-list()
for(LFCthresh in LFCthresholds){
  updown<-NULL
  for (c in contrastNames){
    #sample<-gsub("\\.deseq2\\.results\\.tsv","",t)
    df<-read.delim(paste0(workDir,"tables/differential/",c,".deseq2.results.tsv"),header=T,stringsAsFactors=F)
    df$padj[is.na(df$padj)]<-0 # set NA to 1
    tmp<-data.frame(sample=c,
                     up=sum(df$padj<0.05 & df$log2FoldChange>LFCthresh),
                     down=sum(df$padj<0.05 & df$log2FoldChange<(-LFCthresh)))
    if(is.null(updown)){
      updown<-tmp
    } else {
      updown<-rbind(updown,tmp)
    }
  }
  tbl[[as.character(LFCthresh)]]<-updown
}

sink("./custom/txt/summaryUpDownDiffThresholds.txt")
for (t in 1:length(tbl)) {
  print(paste("padj <0.05 & LFC Threshold:", names(tbl)[t]))
  print(tbl[[t]])
  cat("\n")  # Adds spacing between tables
}
sink()

## interactive Volcano plots -------
results<-readRDS(paste0("custom/rds/",prefix,".results_annotated.RDS"))
dir.create("./custom/plots/volcano", showWarnings = FALSE, recursive = TRUE)
rockman<-readRDS("/Volumes/external.data/MeisterLab/publicData/Various/chrRegions_Rockman2009.RDS")
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
  tmp<-res[res$group==contrasts$shortId[i],c("padj","log2FoldChange","gene_name","gene_biotype","chrRegion","chrRegionType")]
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
  ggsave(filename = paste0("./custom/plots/volcano/volcano_",contrasts$id[i],"_free.png"),
         plot = volcano_plot1, width = 8, height = 6, dpi = 300,bg="white")
  volcano_plot
  interactive_volcano <- ggplotly(volcano_plot,tooltip="text")
  saveWidget(interactive_volcano, file = paste0("./custom/plots/volcano/volcano_",contrasts$id[i],"_free.html"))
}

# zoom LFC -5 to 5
for(i in 1:length(contrasts$id)){
  tmp<-res[res$group==contrasts$shortId[i],c("padj","log2FoldChange","gene_name","gene_biotype","chrRegion","chrRegionType")]
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
  ggsave(filename = paste0("./custom/plots/volcano/volcano_",contrasts$id[i],"_-5to5.png"),
         plot = volcano_plot1, width = 8, height = 6, dpi = 300,bg="white")
  volcano_plot
  interactive_volcano <- ggplotly(volcano_plot,tooltip="text")
  saveWidget(interactive_volcano, file = paste0("./custom/plots/volcano/volcano_",contrasts$id[i],"_-5to5.html"))
}


## Up down by chromosome  -----
dir.create("./custom/plots/byChromosome", showWarnings = FALSE, recursive = TRUE)
nuclear<-seqlevels(Celegans)[1:6]
resA<-res[res$seqnames %in% nuclear,]
resA$padj[is.na(resA$padj)]<-1
padjVal=0.05
lfcVal=0.5
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
ggsave(filename = paste0("./custom/plots/byChromosome/lfcByChr_",
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
ggsave(filename = paste0("./custom/plots/byChromosome/countsByChr_",
                         contrasts$id[i],"_padj",padjVal,"_lfc",lfcVal,".png"),
       plot = p, width = 11, height = 8, dpi = 300,bg="white")




## Up down chromosome regions -----
dir.create("./custom/plots/byChrRegion", showWarnings = FALSE, recursive = TRUE)
autosomes<-seqlevels(Celegans)[1:5]
resA<-res[res$seqnames %in% autosomes,]
resA$padj[is.na(resA$padj)]<-1
resA$chrRegionType<-factor(resA$chrRegionType, levels=c("tip","arm","center"))
padjVal=0.05
lfcVal=0
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
ggsave(filename = paste0("./custom/plots/byChrRegion/lfcByChrRegion_",
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

ggsave(filename = paste0("./custom/plots/byChrRegion/countsByChrRegion_",
                         contrasts$id[i],"_padj",padjVal,"_lfc",lfcVal,".png"),
       plot = p, width = 10, height = 8, dpi = 300,bg="white")


## Heatmap of significant genes -----
results<-readRDS(paste0("./custom/rds/",prefix,".results_annotated.RDS"))
dir.create("./custom/plots/heatmaps", showWarnings = FALSE, recursive = TRUE)
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


getLFCmat<-function(results, numSamplesSignificant=1){
  mat_padj<-gatherResults(results,valueColumn="padj")
  mat_padj[is.na(mat_padj)]<-1
  sigGenes<-mat_padj$gene_id[rowSums(mat_padj[,2:ncol(mat_padj)]<0.05)>numSamplesSignificant]
  mat_lfc<-gatherResults(results,valueColumn="log2FoldChange")
  mat_lfc<-mat_lfc[mat_lfc$gene_id %in% sigGenes,]
  rownames(mat_lfc)<-mat_lfc$gene_id
  mat_lfc<-mat_lfc[,2:ncol(mat_lfc)]
  mat_lfc<-as.matrix(mat_lfc)
  colnames(mat_lfc)<-gsub("log2FoldChange_","",colnames(mat_lfc))
  return(mat_lfc)
}


numSamplesSignificant=3
mat_lfc<-getLFCmat(res_n2,numSamplesSignificant=numSamplesSignificant)

png(paste0("./custom/plots/heatmaps/hclust_heatmap_sigSamples",numSamplesSignificant,"_n2Contrasts.png"),width=19,height=29,units="cm",res=150)
ht<-Heatmap(mat_lfc,show_row_names=F, cluster_columns=T, cluster_rows=T, show_row_dend = F,
            column_names_rot=75,column_title=paste0("log2FC of ",nrow(mat_lfc)," genes singificant in >=",numSamplesSignificant," samples"))
ht<-draw(ht)
dev.off()


mat_lfc<-getLFCmat(res_lin61,numSamplesSignificant=numSamplesSignificant)

png(paste0("./custom/plots/heatmaps/hclust_heatmap_sigSamples",numSamplesSignificant,"_lin61Contrasts.png"),width=19,height=29,units="cm",res=150)
ht<-Heatmap(mat_lfc,show_row_names=F, cluster_columns=T, cluster_rows=T, show_row_dend = F,
            column_names_rot=75,column_title=paste0("log2FC of ",nrow(mat_lfc)," genes singificant in >=",numSamplesSignificant," samples"))
ht<-draw(ht)
dev.off()



numSamplesSignificant=1
mat_lfc<-getLFCmat(results,numSamplesSignificant=numSamplesSignificant)

png(paste0("./custom/plots/heatmaps/hclust_heatmap_sigSamples",numSamplesSignificant,"_allContrasts.png"),width=19,height=29,units="cm",res=150)
ht<-Heatmap(mat_lfc,show_row_names=F, cluster_columns=T, cluster_rows=T, show_row_dend = F,
            column_names_rot=90,column_title=paste0("log2FC of ",nrow(mat_lfc)," genes singificant in >=",numSamplesSignificant," samples"))
ht<-draw(ht)
dev.off()
