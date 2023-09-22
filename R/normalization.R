
#' rewrite scale function, remove all other attributes except normalized data
#'
#' @param input_numeric input_data
#' @param scale_args  args pass to scale function
#'
#' @return normalized data
#' @export
#'
#' @examples
fun_norm_scale <- function(input_numeric,scale_args = list()){
  dims <- dim(input_numeric)
  scale_args$x <- input_numeric
  input_numeric <- do.call("scale",scale_args)
  attr(input_numeric,"scaled:center") <- NULL
  attr(input_numeric,"scaled:scale") <- NULL
  dim(input_numeric) <- dims
  return(input_numeric)
}

#' min max normalization, change data int 0-1
#'
#' @param input_numeric input_numeric data,DO NOT INPUT MATRIX
#'
#' @return normalized numeric
#' @export
#'
#' @examples
fun_norm_minmax <- function(input_numeric){
  input_numeric = input_numeric - min(input_numeric,na.rm = T)
  input_numeric = input_numeric/max(input_numeric,na.rm = T)
  input_numeric
}

#' normalized(been subtracted) input_mat1 by rowMeans of input_mat2
#'
#' @param input_mat1 data used to normalize(kept)
#' @param input_mat2 data used for normalizing(dropped)
#'
#' @return normalized input_mat1
#' @export
#'
#' @examples
fun_norm_substract <- function(input_mat1,input_mat2){
  # 1) check shape
  stopifnot("The nrow of input_mat1 shoud be same with input_mat2" = (all(nrow(input_mat1) == nrow(input_mat2))))
  stopifnot("found emtpy input_mat1" = (ncol(input_mat1) <= 0 ))
  stopifnot("found emtpy input_mat2" = (ncol(input_mat2) <= 0 ))

  # 2) get rowsum of input_mat2
  input_mat2_row_means = as.numeric(rowMeans(input_mat2,na.rm = T))
  apply(input_mat1,2,function(x) x - input_mat2_row_means)

}


#' data.frame normlization
#'
#' @param input_df data.frame obj
#' @param input_variables variables to normalized
#' @param output_variable output_variables for input_varaibles
#' @param method normalization methods, "scale","min-max","towards"
#' @param towards_sample towards_sample,currently not suppored
#' @param fun_norm_scale_args fun_norm_scale_args
#' @param fun_norm_minmax_args fun_norm_minmax_args
#' @param fun_norm_substract_args fun_norm_substract_args
#'
#' @return normlized dataframe
#' @export
#'
#' @examples
fun_norm_df <- function(input_df,
                        input_variables,
                        output_variable=NULL,
                        method =c("scale","min-max","towards"),
                        towards_sample = NULL,
                        fun_norm_scale_args = list(),
                        fun_norm_minmax_args = list(),
                        fun_norm_substract_args = list()){
  input_variables = unique(input_variables)
  stopifnot("some input_variables not found in input_df" = (all(input_variables %in% colnames(input_df))))

  # 1) get methods
  method = match.arg(method,c("scale","min-max","towards"))
  if(method == "scale"){
    current_method_choose_for_normalization = "fun_norm_scale"
    args = fun_norm_scale_args

  }else if(method == "min-max"){
    current_method_choose_for_normalization = "fun_norm_minmax"
    args = fun_norm_minmax_args

  }else if(method == "towards"){
    stop("currently not supported for data.frame")
    current_method_choose_for_normalization = "fun_norm_substract"
    args = fun_norm_substract_args

  }else{
    stop("invalid method")
  }

  # 2) get output_variables
  output_variable =fun_utils_broadcast(
    input_variables = input_variables,
    input_targets = output_variable
  )

  # 3) for loop

  for(i in seq_along(input_variables)){
    each_var = input_variables[i]
    each_var_output = output_variable[i]

    args$input_numeric = input_df[[each_var]]
    input_df[[output_variable]] = do.call(current_method_choose_for_normalization,args = args)
  }
  input_df
}
