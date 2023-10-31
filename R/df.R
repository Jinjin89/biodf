
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
fun_df_feature_map <- function(input_df,input_variables,from,to,
                               new_name = NULL,not_found_replace = NA_character_){
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
  data.frame(old = old_data,new = new_data) %>%
    dplyr::count(old,new) %>%
    print
  input_df
}

#' combine data.frame
#'
#' @param input_df data.frame
#' @param input_variables input variables to collapse
#' @param new_var new var name
#' @param anno_mark the mark separate the variable name and value
#' @param sep the mark separate the variables
#'
#' @return data.frame
#' @export
#'
#' @examples
fun_df_combine_columns <- function(input_df,input_variables,
                                   new_var = "others",
                                   anno_mark = ":",sep = ";"){
  if(!new_var %in% colnames(input_df)){
    input_df[[new_var]] = ""
  }
  for(each_var in input_variables){
    data_tmp = paste0(each_var,anno_mark,input_df[[each_var]])
    data_tmp = paste0(data_tmp,sep)
    input_df[[new_var]] = paste0(input_df[[new_var]],data_tmp)
  }
  input_df[[new_var]]  = stringr::str_remove(input_df[[new_var]],paste0(sep,"$"))
  input_df
}
