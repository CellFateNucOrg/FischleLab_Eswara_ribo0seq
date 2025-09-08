## add annotation to results tables with gene names and locations.
## Combine all results tables into a single .rds object
## extract table of Number of significant up/down regulated genes by different thresholds

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
prefix="ribo0_canonical_geneset_all"

setwd(workDir)
dir.create(paste0(workDir,runName,"/custom/rds"), showWarnings = FALSE, recursive = TRUE)
dir.create(paste0(workDir,runName,"/custom/txt"), showWarnings = FALSE, recursive = TRUE)
genomeVer<-"WS295"

### functions ------
annotate_results <- function(gtf, pattern="\\.deseq2\\.results\\.tsv", path=".", prefix=""){
  fileList<-list.files(paste0(path,"/tables/differential"),pattern=pattern, full.names=T)
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
    write.table(df, paste0(path,"/custom/txt/",file_name, ".deseq2.results_annotated.tsv"),
                row.names=F, quote=F, sep="\t", col.names=T)
  }
}


combine_results <- function(pattern="\\.deseq2\\.results_annotated\\.tsv",path="."){
  fileList <- list.files(paste0(path,"/custom/txt/"), pattern=pattern, full.names=T)
  results <- do.call(rbind, lapply(fileList, read.delim, header=T))
  return(results)
}

#-----------------


gtf<-import(paste0(serverPath,"/publicData/genomes/",genomeVer,"/c_elegans.PRJNA13758.",genomeVer,".canonical_geneset.gtf"))
gtf <- gtf[gtf$type == "gene"]
mcols(gtf)<-mcols(gtf)[c("source","type","gene_id","gene_biotype","gene_name")]
gtf$source<-genomeVer
gtf<-sort(gtf)



## annotate results tables -----

contrasts<-read.csv(paste0(workDir,"/contrasts.csv"),sep=",",header=T)

annotate_results(gtf,path=paste0(workDir,runName),pattern="\\.deseq2\\.results\\.tsv")

results<-combine_results(path=paste0(workDir,runName),pattern="\\.deseq2\\.results_annotated\\.tsv")
results$seqnames<- factor(results$seqnames, levels=c("I", "II", "III", "IV", "V", "X", "MtDNA"))
results$shortId<-results$group
contrasts$shortId<- gsub("_","",contrasts$id)
write.csv(contrasts,paste0(workDir,"/contrasts.csv"),row.names=F,quote=F)

#idx<-match(results$group,contrasts$shortId)
#results$contrast<-contrasts$id[idx]
results$shortId<-factor(results$group, levels=contrasts$shortId)
results$group<-factor(results$group, levels=contrasts$shortId, labels=contrasts$id)

if (file.exists(paste0(workDir,runName,"/custom/rds/",prefix,".results_annotated.RDS"))) {
  response <- readline(prompt = paste(prefix, "file already exists. Overwrite? [y/n]: "))
  if (tolower(response) == "y") {
    saveRDS(results, paste0(workDir,runName,"/custom/rds/",prefix,".results_annotated.RDS"))
  } else {
    print("file not overwritten")
  }
} else {
  saveRDS(results, paste0(workDir,runName,"/custom/rds/",prefix,".results_annotated.RDS"))
}



## Get counts of up/down by different thresholds-----

LFCthresholds<-c(0,0.5,1,1.5,2)

contrastNames<-contrasts$id
tbl<-list()
for(LFCthresh in LFCthresholds){
  updown<-NULL
  for (c in contrastNames){
    #sample<-gsub("\\.deseq2\\.results\\.tsv","",t)
    df<-read.delim(paste0(workDir,runName,"/tables/differential/",c,".deseq2.results.tsv"),header=T,stringsAsFactors=F)
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

sink(paste0(workDir,runName,"/custom/txt/summaryUpDownDiffThresholds.txt"))
for (t in 1:length(tbl)) {
  print(paste("padj <0.05 & LFC Threshold:", names(tbl)[t]))
  print(tbl[[t]])
  cat("\n")  # Adds spacing between tables
}
sink()

