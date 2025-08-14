library(dplyr)
library(ggplot2)



theme_set(
  theme_classic()+
    theme(panel.grid.major = element_blank(),
          panel.grid.minor = element_blank(),
          axis.title.y=ggtext::element_markdown(size=9),
          axis.title.x=ggtext::element_markdown(size=9),
          title=ggtext::element_markdown(size=9)
    )
)

# compare old and new libraries

oldDir<-"/Volumes/external.data/MeisterLab/FischleLab_KarthikEswara/ribo0seq/diff_abund_1_canonical_geneset"
newDir<-"/Volumes/external.data/MeisterLab/FischleLab_KarthikEswara/ribo0seq/diff_abund_2_canonical_noRRnoSP"

oldsf<- read.csv(paste0(oldDir,"/other/deseq2/hpl2_lin61_vs_N2.deseq2.sizefactors.tsv"),
                 stringsAsFactors = FALSE, sep="\t")
newsf<-read.csv(paste0(newDir,"/other/deseq2/hpl2_lin61_vs_N2.deseq2.sizefactors.tsv"),
                stringsAsFactors = F, sep="\t")

sizeFactors<-left_join(oldsf,newsf,by="sample")
sizeFactors$group<-gsub("_Rep.","",sizeFactors$sample)

sizeFactors

ggplot(sizeFactors,aes(x=sizeFactor.x,y=sizeFactor.y,color=group)) +
  geom_point(size=3, alpha=0.7) +
  geom_abline(slope=1, intercept=0, color="black") +
  labs(x="Old size factor", y="New size factor",
       title="Comparison of size factors +- extreme outliers")
