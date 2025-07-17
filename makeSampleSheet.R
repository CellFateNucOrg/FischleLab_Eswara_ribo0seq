workDir<-"/mnt/external.data/MeisterLab/FischleLab_KarthikEswara/ribo0seq/"
ff<-read.delim("./fastqFiles.txt",header=F)
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

#write.csv(df,"./samplesheet_extended.csv",row.names=F,
#          quote=F)
#correctCols<-c("sample","fastq_1","fastq_2","strandedness")


write.csv(df[,correctCols],"./samplesheet.csv",row.names=F,
          quote=F)




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


