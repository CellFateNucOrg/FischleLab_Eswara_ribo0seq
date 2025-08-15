
#' Convert DESeq2 results table to GRanges object
#'
#' @param results DESeq2 results table.
#' @param sort Logical, whether to sort the GRanges object by genomic coordinates. (default: TRUE)
#' @return GRanges object with the same data as the input results table.
#' @export
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

## for heatmaps ------

#' Gather DESeq2 results by different contrasts (groups)
#'
#' Takes a DESeq2 results table joined froom many contrasts with the contrast
#' name in one of the columns and gathers the results into a wide format
#' @param results DESeq2 results table with a column for the contrast name.
#' @param valueColumn Name of the column containing the values to gather (e.g. padj or log2FoldChange)
#' @param nameColumn Name of the column containing the contrast names that will be used as the new column names (default: "group")
#' @return a wide dataframe with gene_id column and value columns.
#' @export
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


#' Make a matrix of particular values from a DESeq2 results table with multiple contrasts (groups)
#'
#' Takes a DESeq2 results table joined froom many contrasts with the contrast
#' name in one of the columns and produces a matrix if LFC values for all genes
#' that are significant in at least numSamplesSignificant samples with a particular
#' log2FoldChange thresholds (lfcVal) and adjusted p value (padjVal)
#' @param results DESeq2 results table with a column for the contrast name.
#' @param numSamplesSignificant Number of samples in which a genes has to be significantly
#' changing in order to be included in the matrix (default: 1).
#' @param lfcVal Log2 fold change threshold for significance (default: 0).
#' @param padjVal Adjusted p-value threshold for significance (default: 0.05).
#' @param direction Direction of regulation to consider: "both", "up", or "down".
#' @param chromosomes Chromosomes to be considered in the analysis: "all" or "autosomal".
#' @return A matrix of LFC values for genes that are significant in at least numSamplesSignificant samples.
#' @export
getLFCMat<-function(results, numSamplesSignificant=1, lfcVal=0, padjVal=0.05,
                    direction="both", chromosomes="all"){
  if(chromosomes=="all"){
    chr<-c("I","II","III","IV","V","X")
    res<-results[results$seqnames %in% chr,]
  } else if(chromosomes=="autosomal"){
    chr<-c("I","II","III","IV","V")
    res<-results[results$seqnames %in% chr,]
  } else {
    stop("Invalid chromosomes specified. Use 'all' or 'autosomal'.")
  }
  mat_padj<-gatherResults(res,valueColumn="padj")
  mat_padj[is.na(mat_padj)]<-1
  mat_lfc<-gatherResults(res,valueColumn="log2FoldChange")
  if(direction=="both"){
    sig<-(mat_padj[,2:ncol(mat_padj)]<padjVal & abs(mat_lfc[,2:ncol(mat_lfc)])>lfcVal)
  } else if(direction=="up"){
    sig<-(mat_padj[,2:ncol(mat_padj)]<padjVal & mat_lfc[,2:ncol(mat_lfc)]>lfcVal)
  } else if(direction=="down"){
    sig<-(mat_padj[,2:ncol(mat_padj)]<padjVal & mat_lfc[,2:ncol(mat_lfc)]<(-lfcVal))
  } else {
    stop("Invalid direction specified. Use 'both', 'up', or 'down'.")
  }

  sigGenes<-mat_padj$gene_id[rowSums(sig)>numSamplesSignificant]
  mat_lfc<-mat_lfc[mat_lfc$gene_id %in% sigGenes,]
  rownames(mat_lfc)<-mat_lfc$gene_id
  mat_lfc<-mat_lfc[,2:ncol(mat_lfc)]
  mat_lfc<-as.matrix(mat_lfc)
  colnames(mat_lfc)<-gsub("log2FoldChange_","",colnames(mat_lfc))
  return(mat_lfc)
}


#' Make heatmap of LFC values
#'
#' Make an heatmapt plot for genes in `res` that are significant in at least `numSamplesSignificant`
#' groups, using a log2 fold change threshold of `lfcVal` and an adjusted p-value threshold of `padjVal`.
#' @param res Results data frame containing gene expression data from multiple samples.
#' @param numSamplesSignificant Minimum number of samples in which a gene must be significant.
#' @param lfcVal Minimum absolute log2 fold change value for significance.
#' @param padjVal Maximum adjusted p-value for significance.
#' @param direction Direction of regulation to consider: "both", "up", or "down".
#' @param groupsToPlot Names of the groups to be included in the upset plot. If NULL, all groups are included.
#' @param setName Name of the set to be included in the plot filename.
#' @param chromosomes Chromosomes to be considered in the analysis: "all" or "autosomal".
#' @return A ggplot object representing the heatmap plot. it is also saved as a PNG file in the ./custom/plots/heatmaps/.
#' @export
makeHeatmapPlot<-function(res, numSamplesSignificant, lfcVal, padjVal, direction,
                        groupsToPlot=NULL, setName="",
                        chromosomes="all"){
  mat_lfc<-getLFCMat(res, numSamplesSignificant, lfcVal, padjVal, direction, chromosomes)
  if(is.null(groupsToPlot)){
    groupsToPlot<-colnames(mat_lfc)
  }
  if(!(setName=="")){
    setName=paste0(setName,"_")
  }
  png(paste0(workDir, runName,"/custom/plots/heatmaps/",setName,"hclust_heatmap_sigSamples",
             numSamplesSignificant,"_padj",padjVal,"_lfc",lfcVal,"_",direction,"_",chromosomes,"Chr.png"),
      width=19,height=29,units="cm",res=150)
  plotTitle<-paste0("Log2FC for ",nrow(mat_lfc)," ",setName," ",chromosomes," chr genes sig. in >=",
                            numSamplesSignificant," samples (padj<",padjVal," ",
                    ifelse(direction=="both","|LFC|","LFC"),
                    ifelse(direction=="down","< -",">"),lfcVal,")")
  ht<-Heatmap(mat_lfc,show_row_names=F, cluster_columns=T, cluster_rows=T, show_row_dend = F,
              column_names_rot=75, column_title=plotTitle,
              col=circlize::colorRamp2(breaks = c(-2, 0, 2), colors = c("blue", "white", "red")))
  ht<-draw(ht)
  dev.off()
  return(ht)
}
