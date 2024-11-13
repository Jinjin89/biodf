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



#'enrichment analysis for input_dta using enriched,
#' @param input_dta genes vector for ora, data.frame for gsea
#' @param input_methods ora or gsea
#' @param outfile the outfile position,default is NULL, not saving the results
#' @param pval pval cutoff for pvalue and qvalue
#' @param logFC logFC
#' @param gene gene
#' @param term2gene term2gene
#' @param return_list return_list
#' @param enrich_db db list, must be names list
#'
#' @return list of different enrichment obj
#' @export
#'

fun_enrich_patch = function(input_dta,input_methods = "ora",
                            outfile = NULL,
                            pval = 0.2,
                            logFC = "logFC",gene = "ID",
                            term2gene = c("term","gene"),
                            enrich_db = list(
                              hallmark = biodata::db_hallmark,
                              kegg = biodata::db_kegg_v2,
                              go = biodata::db_go,
                              reactome = biodata::db_reactome,
                              reactome_metabolism = biodata::db_reactome_metabolism,
                              biocarta = biodata::db_biocarta
                            ),
                            return_list = list()){
  stopifnot('outfile should be null or character' = (is.null(outfile) || is.character(outfile)))
  db_names = names(enrich_db) %>% unique
  db_names =db_names[db_names!='']
  stopifnot('enrich db should have all names' = (length(db_names) == length(enrich_db)))

  if(is.null(outfile) || (!file.exists(outfile))){
    message("enrich results not found run it!")
    enrich_method <- input_methods[1]
    if(enrich_method == 'ora'){
      return_list$degs = data.frame(genes = input_dta)
    }else{
      return_list$degs = input_dta
    }

    for(each_db_name in db_names){
      message("===================================")
      message('enrich with ', enrich_method, ' for ',each_db_name)
      db_now <- enrich_db[[each_db_name]]
      if(tolower(enrich_method) == "ora"){
        return_list[[each_db_name]] <- fun_enrich_ora(input_dta,input_gs = db_now,pval = pval,term2gene = term2gene)
      }else{
        return_list[[each_db_name]] <- fun_enrich_gsea(input_dta,input_gs = db_now,pval = pval,logFC = logFC,gene = gene,term2gene = term2gene)
      }
    }

    # annotate the db
    message("===================================")
    message("Currently only support annotaion for: go, kegg and reactome_metabolism.")
    for(each_db_name in db_names){
      if(each_db_name == 'go'){
        message("annotation GO database with ontology")
        anno_go <-
          biodata::db_go %>% dplyr::count(term,ontology) %>%
          as.data.frame() %>%
          magrittr::set_rownames(.$term)
        return_list$go %<>%
          clusterProfiler::mutate(ontology = anno_go[ID,'ontology'])
      }else if(each_db_name == 'kegg'){
        message("annotation kegg database with class1 and class2")
        anno_kegg <-
          biodata::db_kegg_v2 %>% dplyr::count(term,class1,class2) %>%
          as.data.frame() %>%
          magrittr::set_rownames(.$term)
        return_list$kegg %<>%
          clusterProfiler::mutate(
            class1 = anno_kegg[ID,'class1'],
            class2 = anno_kegg[ID,'class2'])
      }else if(each_db_name == 'reactome_metabolism'){
        message("annotation reactome_metabolism database with pathway")
        anno_reactome_metabolism <-
          biodata::db_reactome_metabolism %>%
          dplyr::count(term,pathway) %>%
          as.data.frame() %>%
          magrittr::set_rownames(.$term)
        return_list$reactome_metabolism %<>%
          clusterProfiler::mutate(pathway = anno_reactome_metabolism[ID,'pathway'])
      }else{
        message("No annotation for: ",each_db_name)
      }
    }

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
