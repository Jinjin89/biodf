#' ORA enerichment using clusterprofiler
#'
#' @param input_genes genes for test
#' @param input_gs geneset data, must have term and gene column
#' @param pval pval cutoff for pvalue and qvalue
#' @param term2gene columns specify term and gene location
#'
#' @return enrichment results
#' @export
#'
#' @examples
fun_enrich_ora = function(input_genes,input_gs,pval = 0.2,term2gene = c("term","gene")){
  # 1) term to gene
  term_2_genes = dplyr::select(input_gs,dplyr::all_of(term2gene))

  # 2) enrich step
  enrich_results = clusterProfiler::enricher(
    gene = input_genes,
    TERM2GENE = term_2_genes,
    pvalueCutoff = pval,
    qvalueCutoff = pval )

  return(enrich_results)

}

#' GSEA analysis
#'
#' @param input_genes data.frame contains gene
#' @param input_gs term2gene
#' @param pval pval cutoff
#' @param logFC logFC column
#' @param gene gene colum
#' @param term2gene columns specify term and gene location
#'
#' @return enrich results
#' @export
#'
#' @examples
fun_enrich_gsea = function(input_genes,input_gs,pval = 0.2,
                           logFC = "logFC",gene = "ID",
                           term2gene = c("term","gene")){
  # 1) term to gene
  term_2_genes = dplyr::select(input_gs,dplyr::all_of(term2gene))


  # 2) prepare the genes info
  input_genes = dplyr::arrange(input_genes,dplyr::desc(!!as.name(logFC)))
  genes = input_genes[[logFC]]
  names(genes) = input_genes[[gene]]

  # 2) enrich step
  enrich_results = clusterProfiler::GSEA(
    geneList = genes,
    TERM2GENE = term_2_genes,
    pvalueCutoff = pval,
    eps=0)

  return(enrich_results)

}



#' GO, kegg and hallmark enrichment analysis for input_data,
#'
#' @param input_df data.frame for gsea or gene list for ora
#' @param input_methods ora or gse
#' @param outfile the outfile postion,default is NULL, not saving the results
#' @param pval pval cutoff for pvalue and qvalue
#' @param logFC
#' @param gene
#' @param term2gene
#' @param return_list
#'
#' @return list of different enrichment obj
#' @export
#'
#' @examples
#' # there is no example
fun_enrich_patch = function(input_df,input_methods = "ora",
                            outfile = NULL,
                            pval = 0.2,
                            logFC = "logFC",gene = "ID",
                            term2gene = c("term","gene"),
                            return_list = list()){
  stopifnot('outfile should be null or character' = (is.null(outfile) || is.character(outfile)))

  if(is.null(outfile) || (!file.exists(outfile))){
    message("enrich results not found run it!")
    go_db = biodata::db_go
    kegg_db = biodata::db_kegg
    hallmark_db = biodata::db_hallmark
    if(tolower(input_methods[1]) == "ora"){
      df_return = data.frame(genes = ora)
      go_res = fun_enrich_ora(input_df,input_gs = go_db,pval = pval,term2gene = term2gene)
      kegg_res = fun_enrich_ora(input_df,input_gs = kegg_db,pval = pval,term2gene = term2gene)
      hallmark_res = fun_enrich_ora(input_df,input_gs = hallmark_db,pval = pval,term2gene = term2gene)
    }else{
      df_return = input_df
      go_res = fun_enrich_gsea(input_df,input_gs = go_db,pval = pval,logFC = logFC,gene = gene,term2gene = term2gene)
      kegg_res = fun_enrich_gsea(input_df,input_gs = kegg_db,pval = pval,logFC = logFC,gene = gene,term2gene = term2gene)
      hallmark_res = fun_enrich_gsea(input_df,input_gs = hallmark_db,pval = pval,logFC = logFC,gene = gene,term2gene = term2gene)
    }
    return_list$degs = df_return
    return_list$go = go_res
    return_list$kegg = kegg_res
    return_list$hallmark = hallmark_res
    if(is.null(outfile)){
      message('Not saving the results into disk!')
    }else{
      message("Saveing results into disk!")
      return_list %>% saveRDS(outfile)
      return_list <- readRDS(outfile)
    }
  }else if(file.exists(outfile)){
    message("enrich results found, load it!")
    return_list = readRDS(outfile)
  }else{
    stop("Something happned, check the source code")
  }
  invisible(return_list)
}
