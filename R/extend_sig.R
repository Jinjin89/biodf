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

#' Merge gene mutation status into clinical data
#'
#' @param input_df imvigor clinical data
#' @param input_variables genes
#'
#' @return data.frame, the merged df with genes mutation status
#' @export
#'
#' @examples
fun_extend_imvigor_mut <- function(
    input_df,input_variables){
  # 1)
  require(IMvigor210CoreBiologies)
  data("fmone")

  # 2) get_genes_mutation-status
  known_short = fmone@assayData$known_short
  known_short = as.data.frame(known_short)
  known_short = as.data.frame(t(known_short))
  known_short$accession = rownames(known_short)
  known_short = tidyr::pivot_longer(known_short,-accession)

  # 3) get samples_patient
  pheno_data = fmone@phenoData@data
  known_short$ID = pheno_data[known_short$accession,"ANONPT_ID"]
  known_short$sample = rownames(input_df)[match(known_short$ID,input_df$ANONPT_ID)]
  #known_short = as.data.frame(known_short)

  known_short$value = ifelse(known_short$value == "",0,1)
  known_short %<>%
    dplyr::filter(!is.na(.$sample)) %>%
    dplyr::filter(value != 0) %>%
    dplyr::count(name,sample)
  all_samples_mut = known_short$sample %>% unique

  if(is.null(input_variables)){
    input_variables = unique(known_short$name)
  }else{
    input_variables = intersect(unique(known_short$name),input_variables)
  }

  for(g in input_variables){
    samples_found =
      known_short %>%
      dplyr::filter(.$name == g) %>%
      dplyr::pull(sample) %>%
      unique

    input_df[[paste0(g,"_status")]] =
      ifelse(input_df$sample %in% samples_found,"mt",
             ifelse(input_df$sample %in% all_samples_mut,"wt",NA))
  }
  input_df

}

