#' expression filter and/or quantitle normalization by row, if you want to do it by column, trannpose it and then transpose back
#'
#' @param input_mat input_mat
#' @param zero_pct zero_pct,values should be in 0,1
#' @param na_pct na_pct values should be in 0,1
#' @param qn quantitle.normalzation by preprocess core
#' @param qn_args args for quantile normalization
#'
#' @return filtered/normalized matrix
#' @export
#'
#' @examples
fun_mat_filter <- function(input_mat, zero_pct = 0.2, na_pct = NULL, qn = F,qn_args = list(keep.names=TRUE)){
  message(paste0("input shape: ",dim(input_mat)[1],"-",dim(input_mat)[2]))
  if (!is.null(zero_pct)) {
    message("filtering rows with zeros:")
    zero_vaue = min(input_mat,na.rm=T)
    minimun_count = ncol(input_mat) * zero_pct
    keep_index = rowSums(input_mat <= zero_vaue) < minimun_count
    input_mat = input_mat[keep_index, ]
    message(paste0("input shape: ",dim(input_mat)[1],"-",dim(input_mat)[2]))
  }
  if (!is.null(na_pct)) {
    message("fitering rows with NAs:")
    minimun_count = ncol(input_mat) * na_pct
    keep_index = rowSums(is.na(input_mat)) < minimun_count
    input_mat = input_mat[keep_index,]
    message(paste0("input shape: ",dim(input_mat)[1],"-",dim(input_mat)[2]))
  }
  if (qn) {
    message("perform quantile normalizations")
    qn_args$x = input_mat
    input_mat = do.call(preprocessCore::normalize.quantiles,qn_args)
  }
  input_mat
}


#' matrix normalization
#'
#' @param input_mat input_matrix
#' @param method normalization methods, "scale","min-max","towards"
#' @param axis 0 means row, 1 means column
#' @param towards_sample sampe_names used for towards
#' @param fun_norm_scale_args fun_norm_scale_args
#' @param fun_norm_minmax_args fun_norm_minmax_args
#' @param fun_norm_substract_args fun_norm_substract_args
#'
#' @return normalized matrix
#' @export
#'
#' @examples
fun_mat_norm <- function(input_mat,method =c("scale","min-max","towards"),axis=1,
                         towards_sample = NULL,
                         fun_norm_scale_args = list(),
                         fun_norm_minmax_args = list(),
                         fun_norm_substract_args = list()){
  stopifnot("input_mat shold be matrix" == is.matrix(input_mat))
  stopifnot("axis only support 0 and 1" == axis %in% c(0,1))
  # 1) choose the methods
  method = match.arg(method,c("scale","min-max","towards"))

  # 2) choose the method
  if(method == "scale"){
    fun_norm_scale_args$input_numeric = input_mat
    # run by column
    if(axis == 1){
      message("scale by column")
      input_mat = do.call(fun_norm_scale,fun_norm_scale_args)
    }else{
      # run by row
      message("scale by row")
      input_mat = t(input_mat)
      input_mat = do.call(fun_norm_scale,fun_norm_scale_args)
      input_mat = t(input_mat)
    }

  }else if(method == "min-max"){
    message("normalized into range(0,1),if found record with all same values, the results would be all NaN")

    # run by column
    if(axis == 1){
      message("min-max by column")
      for(i in seq_along(input_mat[1,])){
        fun_norm_minmax_args$input_numeric = input_mat[,i]
        input_mat[,i] = do.call(fun_norm_minmax,fun_norm_minmax_args)
      }
    }else{
      # run by row
      message("min-max by row")
      for(i in seq_along(input_mat[,1])){
        fun_norm_minmax_args$input_numeric = input_mat[i,]
        input_mat[i,] = do.call(fun_norm_minmax,fun_norm_minmax_args)
      }
    }

  }else if(method == "towards"){
    stopifnot("please supply which sample to subtract when using fun_norm_towards" = (length(towards_sample) == 0))
    if(axis == 1){
      message("towards by column")
      message(paste0("before towards normalization shape: ",dim(input_mat)[1],":",dim(input_mat)[2]))
      input_mat2 = input_mat[,towards_sample]
      keep_sample = dplyr::setdiff(colnames(input_mat),colnames(input_mat2))
      input_mat1 = input_mat[,keep_sample]
      fun_norm_substract_args$input_mat1 = input_mat1
      fun_norm_substract_args$input_mat2 = input_mat2
      input_mat = do.call(fun_norm_substract,fun_norm_substract_args)
      message(paste0("after towards normalization shape: ",dim(input_mat)[1],":",dim(input_mat)[2]))

    }else{
      message("towards by row")
      message(paste0("before towards normalization shape: ",dim(input_mat)[1],":",dim(input_mat)[2]))
      input_mat = t(input_mat) # transponse: row to column
      input_mat2 = input_mat[,towards_sample]
      keep_sample = dplyr::setdiff(colnames(input_mat),colnames(input_mat2))
      input_mat1 = input_mat[,keep_sample]
      fun_norm_substract_args$input_mat1 = input_mat1
      fun_norm_substract_args$input_mat2 = input_mat2
      input_mat = do.call(fun_norm_substract,fun_norm_substract_args)
      input_mat = t(input_mat) # transponse: column to row(keep the orginal data shape)
      message(paste0("after towards normalization shape: ",dim(input_mat)[1],":",dim(input_mat)[2]))
    }
  }else{
    stop("not found valid input please check input_matrix or method arguments")
  }

  return(input_mat)
}



fun_mat_bind_rows <- function(x,y,keep = c("all","x","y"),x_name = "_x",y_name = "_y",fill_value=NA){
  message("this function could be extremly slow using bigger data, as it use for loop to access each element")
  stopifnot("all types of x and y should be eaqul" = (all(typeof(x) == typeof(y))))
  stopifnot("x and y should be matrix,and should have names,if no name avaiable, use rbind or cbind" = (is.matrix(x) && is.matrix(y) && (!is.null(names(x))) && (!is.null((names(y))))))
  keep = match.arg(keep, c("all","x","y"))

  # 1) get rownames
  x_rowname = rownames(x)
  y_rowname = rownames(y)

  x_colname = colnames(x)
  y_colname = colnames(y)
  common_row_names =  dplyr::intersect(x_rowname,y_rowname) # used for generating new file

  # 3) how to select data
  if(keep == "all"){
    message("keep all the data in x and y")
    if(length(common_row_names) == 0){
      message("no duplicated data found")
      row_names = c(x_rowname,y_rowname)
      col_names = c(x_colname,y_colname)
      # init data
    }else{
      message("duplicated data found")
    }

  }else if(keep =="x"){
    message("keep all the data in the x, append new data in y at the bottom")

  }else if(keep == "y"){

  }else{
    stop("wrong arguments")
  }
}


fun_mat_bind_cols <- function(x,y,keep = c("all","x","y")){

}
