
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


