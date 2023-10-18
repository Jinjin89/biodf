
#' relabel features with new label in element-wise way
#'
#' @param input_df data.frame
#' @param input_variables the variable name in data.frame
#' @param from the old element label
#' @param to the new element label
#' @param new_name new_name, if not supplied, change inplace
#' @param not_found_replace the label
#'
#' @return relabeled dataframe for supplied variable
#' @export
#'
#' @examples
#' fun_df_feature_map(biodata::tcga_clin_xena,
#'   input_variables = "gender",
#'   from = c("MALE","FEMALE"),
#'   to = c("M","F"),
#'   new_name = "sex")
fun_df_feature_map <- function(input_df,input_variables,from,to,new_name = NULL,not_found_replace = NA_character_){
  # 1)
  stopifnot("Input length should be eaual" = (length(from) == length(to)))

  # 2) all the input data summary
  old_data = input_df[[input_variables]]
  input_all_data = unique(input_df[[input_variables]])
  data_not_found_in_from <- dplyr::setdiff(input_all_data,from)
  if(is.null(new_name)){
    message("new name is null, change data inplace")
    new_name = input_variables
  }
  input_df[[new_name]] = input_df[[input_variables]]
  if(length(data_not_found_in_from) > 0){
    warning("Some data not supplied: ", paste0(data_not_found_in_from,", "),
            "set this data with no_found_replace(defualt is NA) ")
    input_df[[new_name]] = ifelse(
      input_df[[new_name]] %in% data_not_found_in_from,
      not_found_replace,
      input_df[[new_name]]
    )
  }
  for(i in seq_along(from)){
    input_df[[new_name]] =
      ifelse(input_df[[new_name]] == from[i],to[i],input_df[[new_name]])
  }
  new_data = input_df[[new_name]]
  print(table(old_data,new_data))
  input_df
}
