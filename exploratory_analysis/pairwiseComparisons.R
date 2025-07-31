library(GenomicRanges)
library(rtracklayer)
library(ggplot2)
library(dplyr)
library(BSgenome.Celegans.UCSC.ce11)
library(ggpubr)
library(plotly)
library(ggrepel)
library(rstatix)


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


## correlations -----
# N2 contrasts
results<-readRDS(paste0(workDir, runName,"/custom/rds/",prefix,".results_annotated.RDS"))
dir.create(paste0(workDir, runName,"/custom/plots/correlation"), showWarnings = FALSE, recursive = TRUE)

n2contrasts<-grep("_vs_N2$",levels(results$group))
lin61contrasts<-grep("_vs_HPL2GFP__lin61$",levels(results$group))

combn(n2contrasts,2)
combn(lin61contrasts,2)
