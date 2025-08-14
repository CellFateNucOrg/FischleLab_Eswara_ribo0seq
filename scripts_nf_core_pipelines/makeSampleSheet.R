library(GenomicRanges)
library(rtracklayer)
library(BSgenome.Celegans.UCSC.ce11)

workDir<-"/Volumes/external.data/MeisterLab/FischleLab_KarthikEswara/ribo0seq/"
dir.create(workDir,showWarnings=F, recursive=T)

ff<-read.delim(paste0(workDir,"/fastqFiles.txt"),header=F)
ff$sample<-sapply(strsplit(basename(ff$V1),"_"),"[[",1)
#ff$replicate<-gsub("Rep","R",sapply(strsplit(basename(ff$V1),"_"),"[[",2))
ff$replicate<-sapply(strsplit(basename(ff$V1),"_"),"[[",2)
ff$R1vR2<-gsub("\\.fq\\.gz","",sapply(strsplit(basename(ff$V1),"_"),"[[",3))
metadata<-data.frame(sample=c("EM88","EM90","EM92","EM91","Btm", "N2","EM38","hpl2tm","lin61tm"),
                     genotype=c("HPL2GFP_lin61","HPL2GFPneutHng_lin61","HPL2GFP3xhngI158A_lin61",
                                "HPL2GFPI158A_lin61","hpl2_lin61", "N2","HPL2GFP", "hpl2","lin61"),
                     phenotype=c("HPL2GFP_lin61","onlyDimer_lin61","onlyLLPS_lin61",
                                 "noDimerLLPS_lin61","hpl2_lin61", "N2","HPL2GFP", "hpl2","lin61"))
metadata
ff$phenotype<-metadata$phenotype[match(ff$sample,metadata$sample)]
ff$genotype<-metadata$genotype[match(ff$sample,metadata$sample)]

df<-data.frame(sample=paste0(ff$phenotype[ff$R1vR2=="1"],"_",ff$replicate[ff$R1vR2=="1"]),
              fastq_1=paste0(workDir,ff$V1[ff$R1vR2=="1"]),
              fastq_2=paste0(workDir,ff$V1[ff$R1vR2=="2"]),
              strandedness="auto",
              sampleID=ff$sample[ff$R1vR2=="1"],
              phenotype=ff$phenotype[ff$R1vR2=="1"],
              genotype=ff$genotype[ff$R1vR2=="1"],
              replicate=ff$replicate[ff$R1vR2=="1"])

write.csv(df,"./samplesheet_extended.csv",row.names=F,
          quote=F)
#correctCols<-c("sample","fastq_1","fastq_2","strandedness")


#write.csv(df[,correctCols],"./samplesheet.csv",row.names=F,
#          quote=F)




contrasts<-data.frame(id=c(paste0(c("HPL2GFP_lin61","onlyDimer_lin61","onlyLLPS_lin61",
                                    "noDimerLLPS_lin61","hpl2_lin61", "HPL2GFP", "hpl2","lin61"),"_vs_N2"),
                           paste0(c("lin61","onlyDimer","onlyLLPS",
                                    "noDimerLLPS","hpl2"),"_vs_HPL2GFP__lin61")),
                      variable="phenotype",
                      reference=c(rep("N2",8),rep("HPL2GFP_lin61",5)),
                      target=c(c("HPL2GFP_lin61","onlyDimer_lin61","onlyLLPS_lin61",
                               "noDimerLLPS_lin61","hpl2_lin61", "HPL2GFP", "hpl2","lin61"),
                               c("lin61","onlyDimer_lin61","onlyLLPS_lin61",
                                 "noDimerLLPS_lin61","hpl2_lin61")),
                      blocking="replicate")

write.table(contrasts,file="./contrasts.csv",sep=",",row.names=F, quote=F)


## remove extreme outliers
genomeVer<-"WS295"
gtf<-import(paste0("/Volumes/MeisterLab/publicData/genomes/",genomeVer,"/c_elegans.PRJNA13758.",genomeVer,".canonical_geneset.gtf"))
#gtf<-import("Z:/MeisterLab/publicData/genomes/WS295/c_elegans.PRJNA13758.WS295.canonical_geneset.gtf")
gtf <- gtf[gtf$type == "gene"]
mcols(gtf)<-mcols(gtf)[c("source","type","gene_id","gene_biotype","gene_name")]
gtf$source<-genomeVer
gtf<-sort(gtf)

# find major ribosomal RNA clusters
idx<-grep("rrn",gtf$gene_name)
rrr<-gtf[idx,]
rrr[seqnames(rrr)=="I"]
rrr[seqnames(rrr)=="V"]




riboRNA_I<-GRanges("I:15060299-15071033")
riboRNA_V<-GRanges("V:17115526-17132107")

rrI<-subsetByOverlaps(gtf,riboRNA_I,ignore.strand=T)
rrI
rrV<-subsetByOverlaps(gtf,riboRNA_V,ignore.strand=T)
rrV
sp<-gtf[grep("srpr",gtf$gene_name)]
sp

toRemove<-c(rrI$gene_id,rrV$gene_id,sp$gene_id)
length(toRemove) # 42 genes

gtf_original<-import(paste0("/Volumes/MeisterLab/publicData/genomes/",genomeVer,"/c_elegans.PRJNA13758.",genomeVer,".canonical_geneset.gtf"))
#gtf_original <- gtf_original[gtf_original$type == "gene"]
norr_nosp<-gtf_original[!(gtf_original$gene_id %in% toRemove),]
export(norr_nosp,paste0("/Volumes/MeisterLab/publicData/genomes/",genomeVer,"/c_elegans.PRJNA13758.",genomeVer,".canonical_geneset_noRR_noSP.gtf"))

# counts<- read.delim(paste0(workDir,"/star_salmon/salmon.merged.gene_counts.tsv"))
# dim(counts)
# counts<-counts[counts$gene_id %in% norr_nosp$gene_id,]
# dim(counts)
# write.table(counts,paste0(workDir,"/star_salmon/salmon.merged.gene_counts_noRR_noSP.tsv"),sep="\t",row.names=F,quote=F)

blacklist<-sort(c(rrI,rrV,sp))
blacklist$name<-blacklist$gene_id
dir.create(paste0(workDir,"/custom/tracks/"),showWarnings=F, recursive=T)
export(blacklist,paste0(workDir,"/custom/tracks/c_elegans.PRJNA13758.",genomeVer,".blacklist_RR_SP.bed"))

seqlevelsStyle(blacklist)<-"UCSC"
blacklist
export(blacklist,paste0(workDir,"/custom/tracks/c_elegans.PRJNA13758.",genomeVer,".blacklist_RR_SP__ce11.bed"))


