library(GenomicRanges)
library(rtracklayer)
library(dplyr)
library(BSgenome.Celegans.UCSC.ce11)

gff2<-import("/Users/semple/Downloads/c_elegans.PRJNA13758.WS295.annotations.gff2.gz")

gff3<-import("/Users/semple/Downloads/c_elegans.PRJNA13758.WS295.annotations.gff3.gz")

genome<-Celegans
dfamVer="Dfam_3.9"

# hmm="https://dfam.org/releases/current/families/Dfam-curated_only-1.hmm.gz"
# ce10 alignments "https://www.dfam.org/releases/annotations/ce10_df38/assembly_alignments/
# https://www.repeatmasker.org/genomes/ce10/RepeatMasker-rm405-db20140131/ce10.fa.out.gz
# https://www.repeatmasker.org/genomes/ce10/RepeatMasker-rm406-dfam2.0/ce10.fa.out.gz
# www.dfam.org/api/families?format=hmm&clade=6239&clade_relatives=both&download=true



args<-commandArgs(trailingOnly=TRUE)
workDir<-args[1]
print(paste0("workDir is: ", workDir))
#workDir="~/Documents/MeisterLab/otherPeopleProjects/Moushumi/SMC_RNAseq_repeats"
workDir<-"/Volumes/external.data/MeisterLab/FischleLab_KarthikEswara/ribo0seq"

dfamURL=paste0("https://www.dfam.org/releases/",dfamVer,"/annotations/ce10/ce10_dfam.nrph.hits.gz")
download.file(url=dfamURL, destfile=paste0(workDir,"/ce10_",dfamVer,".nrph.hits.gz"),method="wget")
system(paste0("gunzip ",workDir,"/ce10_",dfamVer,".nrph.hits.gz"))
repeats_dfam <- read.delim(paste0(workDir,"/ce10_",dfamVer,".nrph.hits"))


repeats.gr <- GRanges(seqnames=repeats_dfam[,"X.seq_name"],
                      ranges=IRanges(start=apply(repeats_dfam[,c("env.st",
                                                                 "env.en")],1,FUN=min),
                                     end=apply(repeats_dfam[,c("env.st",
                                                               "env.en")],1,FUN=max)),
                      strand=repeats_dfam[,"strand"])

repeats.gr$source<-"Dfam_3.3"
repeats.gr$type<-"exon"

repeats.gr$ID<-NA
repeats.gr$Name<-repeats_dfam$family_name
repeats.gr$Alias<-repeats_dfam$family_acc

#export.gff3(repeats.gr, "repeats_ce10_dfam_nr.gff3")

liftoverURL="http://hgdownload.soe.ucsc.edu/goldenPath/ce10/liftOver/ce10ToCe11.over.chain.gz"
download.file(url=liftoverURL,destfile=paste0(workDir,"/ce10ToCe11.over.chain.gz"),method="wget")
system(paste0("gunzip ",workDir,"/ce10ToCe11.over.chain.gz"))

chain <- import(paste0(workDir,"/ce10ToCe11.over.chain"))
# liftover splits some granges, so need to fuse them again
repeats_ce11.gr <- unlist(range(liftOver(repeats.gr,chain)))
table(width(repeats_ce11.gr)-width(repeats.gr)) #checking the differences are minor
mcols(repeats_ce11.gr)<-mcols(repeats.gr)
repeats_ce11.gr$ID<-paste0("rpt",1:length(repeats_ce11.gr),"_",
                           seqnames(repeats_ce11.gr),":",
                           start(repeats_ce11.gr),"-",end(repeats_ce11.gr))

saveRDS(repeats_ce11.gr,paste0(workDir,"/repeats_ce11_",dfamVer,"_nr.rds"))

# bedgraph for Todor
repeats_bg<-repeats_ce11.gr
mcols(repeats_bg)<-NULL
repeats_bg$score<-1
repeats_bg$name<-repeats_ce11.gr$Name
export.bedGraph(repeats_bg,paste0(workDir,"/repeats_ce11_",dfamVer,"_nr.bedgraph"))
export.bed(repeats_bg,paste0(workDir,"/repeats_ce11_",dfamVer,"_nr.bed"))
tileWidth=5000
tiles<-unlist(tileGenome(seqinfo(Celegans),tilewidth=tileWidth))
tiles$score<-countOverlaps(tiles,repeats_bg,ignore.strand=TRUE)
export(tiles,paste0(workDir,"/repeats_ce11_",dfamVer,"_nr_",tileWidth,".bw"),
       format="bigwig")

# make wormbase formatted chr names for gff3
seqlevelsStyle(genome)<-"ensembl"
seqlevelsStyle(repeats_ce11.gr)<-"ensembl"
seqinfo(repeats_ce11.gr)<-seqinfo(genome)
genome(repeats_ce11.gr)<-NA
export.gff3(repeats_ce11.gr, paste0(workDir,"/repeats_ce11_",dfamVer,"_nr.gff3"))

#file.remove(paste0(workDir,"/ce10ToCe11.over.chain.gz"))
file.remove(paste0(workDir,"/ce10ToCe11.over.chain"))
#file.remove(paste0(workDir,"/ce10_dfam.nrph.hits.gz"))
file.remove(paste0(workDir,"/ce10_",dfamVer,".nrph.hits"))
#file.remove(paste0(workDir,"/repeats_ce10_dfam_nr.gff3"))


# get repeat sequences
gff<-import(paste0(workDir,"/repeats_ce11_",dfamVer,"_nr.gff3"),format="gff3")

gtfdf<-data.frame(seqname=seqnames(gff),
                  source=gff$source,
                  feature=gff$type,
                  start=start(gff),
                  end=end(gff),
                  score=".",
                  strand=strand(gff),
                  frame=".",
                  attribute=paste0('gene_id "',gff$ID,'"; gene_name "',
                                   gff$Name,'";'))

write.table(gtfdf,paste0(workDir,"/repeats_ce11_",dfamVer,"_nr.gtf"),sep="\t",
            col.names=F,row.names=F, quote=F)


seqlevelsStyle(genome)<-"ensembl"
rptSeq<-getSeq(genome,gff)

names(rptSeq)<-paste(gff$ID,gff$Name,gff$Alias,sep=";")

writeXStringSet(rptSeq, paste0(workDir,"/repeats_ce11_",dfamVer,"_nr.fa.gz"),
                append=FALSE, compress="gzip", format="fasta")


#### make chrom.sizes
chromSizes<-data.frame(seqnames=seqnames(genome),seqlengths=seqlengths(genome))
write.table(chromSizes,file=paste0(workDir,"/ws235.chrom.sizes"),col.names=F,
            row.names=F,sep="\t",quote=F)



# McMurchy2017 types
fig2data2URL<-"https://elifesciences.org/download/aHR0cHM6Ly9jZG4uZWxpZmVzY2llbmNlcy5vcmcvYXJ0aWNsZXMvMjE2NjYvZWxpZmUtMjE2NjYtZmlnMi1kYXRhMi12My54bHN4/elife-21666-fig2-data2-v3.xlsx?_hash=Vc0gZhYRUq5G8oMoQTJHwwfkDZsgkSAv1k2lEThtEHY%3D"
mcmurchyFilename="elife-21666-fig2-data2-v3.xlsx"

download.file(url=fig2data2URL, destfile=paste0("./",mcmurchyFilename))

mcmurchy<-readxl::read_excel(paste0("./",mcmurchyFilename),sheet="Compiled")

gff<-import(paste0(workDir,"/repeats_ce11_",dfamVer,"_nr.gff3"),format="gff3")

gff$repType<-NA
inMcmurchy<-gff$Name %in% mcmurchy$Family
idxMM<-match(gff$Name[inMcmurchy],mcmurchy$Family)
gff$repType[inMcmurchy]<-mcmurchy$Class[idxMM]

saveRDS(gff,paste0(workDir,"/repeats_ce11_",dfamVer,"_nr.rds"))

write.csv(gff,paste0(workDir,"/repeats_ce11_",dfamVer,"_nr.csv"), quote=F,
          row.names=F)

