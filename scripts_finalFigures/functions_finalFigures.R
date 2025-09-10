
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



#' Get a binary matrix of significant results by group
#'
#' From a results table, extract genes that are signficant in at least
#' `numSamplesSignificant` samples, with an absolute log2 fold change greater
#' than `lfcVal` and adjusted p-value less than `padjVal`. Log2FoldChange values
#' can be evaluated in both directions (up and down), or just up or down
#' regulated.
#' @param results A data frame containing the results of differential expression analysis.
#' @param numSamplesSignificant Minimum number of samples in which a gene must be significant.
#' @param lfcVal Minimum absolute log2 fold change value for significance.
#' @param padjVal Maximum adjusted p-value for significance.
#' @param direction Direction of regulation to consider: "both", "up", or "down".
#' @param chromosomes Chromosomes to be considered in the analysis: "all" or "autosomal".
getBinaryMat<-function(results, numSamplesSignificant=1, lfcVal=0, padjVal=0.05,
                       direction="both", chromosomes="all"){
  if(chromosomes=="all"){
    chr<-c("I","II","III","IV","V","X")
    res<-results[results$seqnames %in% c(chr, paste0("chr",chr)),]
  } else if(chromosomes=="autosomal"){
    chr<-c("I","II","III","IV","V")
    res<-results[results$seqnames %in% c(chr, paste0("chr",chr)),]
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
  rownames(sig)<-mat_padj[,1]
  colnames(sig)<-gsub("padj_","",colnames(sig))
  sig<-sig[rowSums(sig)>=numSamplesSignificant,]

  return(sig)
}
