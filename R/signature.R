#' ssGSEA analysis for input_mat using input_genes_list, the results would be saved as outfile
#'
#' @param input_mat input_mat
#' @param input_genes_list input_genes_list
#' @param outfile outfile
#' @param rds whether save as rds
#'
#' @return data.frame
#' @export
#'
#' @examples

fun_sig_ssGSGA <- function (input_mat, input_genes_list, outfile,rds=F)
{
  if (!file.exists(outfile)) {
    gsea_out = GSVA::gsva(input_mat, input_genes_list, method = "ssgsea",
                          kcdf = "Gaussian", abs.ranking = TRUE, verbose = FALSE)
    gsea_out = data.frame(t(gsea_out),check.names = F)
    gsea_out = dplyr::mutate_all(gsea_out, scale)
    gsea_out = dplyr::mutate_all(gsea_out, as.numeric)
    if(rds){
      saveRDS(gsea_out,outfile)
    }else{
      data.table::fwrite(gsea_out,file = outfile,row.names = T,quote = T)
    }
  }
  if(rds){
    readRDS(outfile)
  }else{
    return(read.csv(outfile,row.names = 1,check.names = F))
  }
}

#' Calculated the weight sum of the given input coefficients
#'
#' @param input_df input_df
#' @param input_coef input_coef
#' @param coef_column which coef to use
#' @param output_name the output name
#'
#' @return data.frame with calculated value
#' @export
#'
#' @examples
fun_sig_weighted_sum = function(
    input_df,
    input_coef,
    coef_column = 1,
    output_name = "Risk_score"){
  # 1) get commom variables
  common_variables = intersect(
    rownames(input_coef),
    colnames(input_df))

  stopifnot("No common variables found"=length(common_variables) > 0)
  message(paste0("Found ",length(common_variables), " variables from ", nrow(input_coef)))
  # 2) getsum
  input_df[[output_name]] = 0
  for(g in common_variables){
    # 1) get current weighted sum
    value_tmp =
      as.numeric(input_df[[g]]) * as.numeric(input_coef[g,coef_column])
    input_df[[output_name]] =
      input_df[[output_name]] + value_tmp
  }

  input_df
}

#' Immune deonvolution
#'
#' @param input_expr expression matrix, should not be log2 transformed
#' @param outfile outfile postion
#' @param indications cancer type
#' @param tumor T or F
#' @param arrays T or F
#'
#' @return data.table
#' @export
#'
fun_sig_deonvolute_immune =  function(input_expr,outfile,indications,tumor=T,arrays=F,perm=1000){

  if(!file.exists(outfile)){
    suppressMessages(require(immunedeconv))
    suppressMessages(require(e1071))
    suppressMessages(require(preprocessCore))
    suppressMessages(require(parallel))
    suppressMessages(require(tidyverse))
    # 1) get celltype map
    cell_map = immunedeconv::cell_type_map
    cell_map = as.data.frame(cell_map)
    cell_map = dplyr::filter(cell_map,
                             method_dataset == "cibersort")
    rownames(cell_map) = cell_map$method_cell_type
    col_names = c("cell_type",colnames(input_expr))

    # 2) TIMER
    if(length(indications) != ncol(input_expr)){
      indications = rep(indications,ncol(input_expr))
    }
    timer_res = immunedeconv::deconvolute(
      gene_expression  = input_expr,
      method = "timer",
      indications = indications
    )
    timer_res$cell_type = paste0(timer_res$cell_type,"_TIMER")

    # 3) CIBERSORT
    message(">>>CIBERSORT")
    lm22 = biodata::sig_list$lm22
    cibersort_res = CIBERSORT(mixture_file = as.data.frame(input_expr),
                              sig_matrix = lm22,
                              perm = perm,
                              QN=arrays,
                              absolute = F)
    cibersort_res =
      cibersort_res %>%
      t %>%
      as.data.frame() %>%
      dplyr::filter(rownames(.) %in% rownames(cell_map)) %>%
      dplyr::mutate(cell_type = cell_map[rownames(.),"cell_type"]) %>%
      dplyr::select(dplyr::all_of(col_names))
    cibersort_res$cell_type = paste0(cibersort_res$cell_type,"_CIBERSORT")


    # 4) cibersort-abs
    message(">>>CIBERSORT-ABS")
    cibersort_abs_res = CIBERSORT(mixture_file = as.data.frame(input_expr),
                                  sig_matrix = lm22,
                                  perm = perm,
                                  QN = arrays,
                                  absolute = T)
    cibersort_abs_res =
      cibersort_abs_res %>%
      t %>%
      as.data.frame() %>%
      dplyr::filter(rownames(.) %in% rownames(cell_map)) %>%
      dplyr::mutate(cell_type = cell_map[rownames(.),"cell_type"]) %>%
      dplyr::select(dplyr::all_of(col_names))
    cibersort_abs_res$cell_type = paste0(cibersort_abs_res$cell_type,"_CIBERSORT-ABS")


    # 5) quantiseq
    quantiseq_res = immunedeconv::deconvolute(
      gene_expression  = input_expr,
      method = "quantiseq",
      arrays = arrays
    )
    quantiseq_res$cell_type = paste0(quantiseq_res$cell_type,"_QUANTISEQ")


    # 6) mcpcounter
    mcpcounter_res = immunedeconv::deconvolute(
      gene_expression  = input_expr,
      method = "mcp_counter"
    )

    mcpcounter_res$cell_type =  paste0(mcpcounter_res$cell_type,"_MCPCOUNTER")

    # 8) xcell
    xcell_res = immunedeconv::deconvolute(
      gene_expression  = input_expr,
      method = "xcell",
      arrays = arrays
    )
    xcell_res$cell_type =  paste0(xcell_res$cell_type,"_XCELL")

    # 9) epic
    epic_res = immunedeconv::deconvolute(
      gene_expression  = input_expr,
      method = "epic",
      tumor = tumor
    )
    epic_res$cell_type =  paste0(epic_res$cell_type,"_EPIC")


    # 10) estimate
    estimate_res = immunedeconv::deconvolute(
      gene_expression  = input_expr,
      method = "estimate"
    )

    estimate_res$cell_type =  paste0(estimate_res$cell_type,"_ESTIMATE")

    estimate_res$cell_type = stringr::str_replace_all(estimate_res$cell_type,"stroma ","Stromal_")
    estimate_res$cell_type = stringr::str_replace_all(estimate_res$cell_type,"immune ","Immune_")
    estimate_res$cell_type = stringr::str_replace_all(estimate_res$cell_type,"estimate ","ESTIMATE_")
    estimate_res = dplyr::filter(estimate_res,cell_type != "tumor purity_ESTIMATE")


    # bindcols
    merge_result =
      dplyr::bind_rows(
        timer_res,
        cibersort_res,
        cibersort_abs_res,
        quantiseq_res,
        mcpcounter_res,
        xcell_res,
        epic_res,
        estimate_res
      )
    merge_result %>%
      as.data.frame() %>%
      set_rownames(.$cell_type) %>%
      dplyr::select(-cell_type) %>%
      t %>%
      as.data.frame() %>%
      dplyr::mutate(sample =rownames(.)) %>%
      select(dplyr::all_of(unique(c("sample",colnames(.))))) %>%
      data.table::fwrite(outfile)
  }
  data.table::fread(outfile)
}
