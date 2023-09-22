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
#' @param pval pval cutoff for pvalue and qvalue
#'
#' @return list of different enrichment obj
#' @export
#'
#' @examples
#'
fun_enrich_patch = function(input_df,input_methods = "ora",pval = 0.2,
                            logFC = "logFC",gene = "ID",
                            term2gene = c("term","gene")){
  go_db = biodata::db_go
  kegg_db = biodata::db_kegg
  hallmark_db = biodata::db_hallmark
  if(tolower(input_methods[1]) == "ora"){
    go_res = fun_enrich_ora(input_df,input_gs = go_db,pval = pval,term2gene = term2gene)
    kegg_res = fun_enrich_ora(input_df,input_gs = kegg_db,pval = pval,term2gene = term2gene)
    hallmark_res = fun_enrich_ora(input_df,input_gs = hallmark_db,pval = pval,term2gene = term2gene)

  }else{
    go_res = fun_enrich_gsea(input_df,input_gs = go_db,pval = pval,logFC = logFC,gene = gene,term2gene = term2gene)
    kegg_res = fun_enrich_gsea(input_df,input_gs = kegg_db,pval = pval,logFC = logFC,gene = gene,term2gene = term2gene)
    hallmark_res = fun_enrich_gsea(input_df,input_gs = hallmark_db,pval = pval,logFC = logFC,gene = gene,term2gene = term2gene)
  }

  return(list(
    go = go_res,
    kegg = kegg_res,
    hallmark = hallmark_res
  ))

}
