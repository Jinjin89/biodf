#' ssGSEA analysis for input_mat using input_genes_list, the results would be saved as outfile
#'
#' @param input_mat input_mat
#' @param input_genes_list input_genes_list
#' @param outfile
#'
#' @return data.frame
#' @export
#'
#' @examples
fun_sig_ssGSGA = function(
    input_mat,
    input_genes_list,
    outfile){
  if(!file.exists(outfile)){
    gsea_out = GSVA::gsva(
      input_mat,
      input_genes_list,
      method="ssgsea",
      kcdf="Gaussian",
      abs.ranking=TRUE,
      verbose=FALSE)
    print(gsea_out)
    gsea_out = data.frame(t(gsea_out))
    gsea_out = dplyr::mutate_all(gsea_out,scale)
    gsea_out = dplyr::mutate_all(gsea_out,as.numeric)
    write.csv(gsea_out,outfile)
  }

  return(read.csv(outfile,row.names = 1))
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
