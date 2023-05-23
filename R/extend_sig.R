#' Calculate the EMT signature
#'
#' @param input_mat matrix
#' @param outfile output file
#'
#' @return data.frame with score
#' @export
#'
#' @examples
fun_extend_sig_emt = function(input_mat,outfile){


  ecm_pos_neg = list()
  ecm_pos_neg$pos = go_db$gene[go_db$term == "GOBP_POSITIVE_REGULATION_OF_EPITHELIAL_TO_MESENCHYMAL_TRANSITION"]
  ecm_pos_neg$neg = go_db$gene[go_db$term == "GOBP_NEGATIVE_REGULATION_OF_EPITHELIAL_TO_MESENCHYMAL_TRANSITION"]

  ecm_df = fun_sig_ssGSGA(
    input_mat =input_mat,
    input_genes_list = ecm_pos_neg,
    outfile = outfile)

  ecm_df$pos46 = ecm_df$pos * 46
  ecm_df$neg28 = ecm_df$neg * 28

  ecm_df$score = ecm_df$pos46 - ecm_df$neg28
  return(ecm_df)

}
