library(GenomicRanges)
library(rtracklayer)
library(BSgenome.Celegans.UCSC.ce11)

workDir<-"/Volumes/external.data/MeisterLab/FischleLab_KarthikEswara/ribo0seq/"
genomeVer<-"WS295"
norr_nosp<-import(paste0("/Volumes/MeisterLab/publicData/genomes/",genomeVer,"/c_elegans.PRJNA13758.",genomeVer,".canonical_geneset_noRR_noSP.gtf"))

counts<- read.delim(paste0(workDir,"/star_salmon/salmon.merged.gene_counts.tsv"))
dim(counts)
counts<-counts[counts$gene_id %in% norr_nosp$gene_id,]
dim(counts)
#write.table(counts,paste0(workDir,"/star_salmon/salmon.merged.gene_counts_noRR_noSP.tsv"),sep="\t",row.names=F,quote=F)
write.table(counts,paste0(workDir,"/star_salmon/salmon.merged.gene_counts_noRR_noSP_moreSeq.tsv"),sep="\t",row.names=F,quote=F)
