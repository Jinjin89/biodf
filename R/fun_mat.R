#' expression filter and/or quantitle normalization by row, if you want to do it by column, trannpose it and then transpose back
#'
#' @param input_mat input_mat
#' @param zero_pct zero_pct[0,1]
#' @param na_pct na_pct [0,1]
#' @param qn quantitle.normalzation by preprocess core
#'
#' @return filtered matrix
#' @export
#'
#' @examples
fun_mat_filter = function(input_mat,zero_pct=0.2,na_pct = NULL,qn = F){
  if(!is.null(zero_pct)){
    # do zero pct
    zero_vaue = min(input_mat)
    minimun_count = ncol(input_mat) * zero_pct

    keep_index = rowSums(input_mat<=zero_vaue) < minimun_count

    input_mat = input_mat[keep_index,]
  }

  if(!is.null(na_pct)){

  }

  if(qn){
    input_mat = preprocessCore::normalize.quantiles(input_mat,keep.names = T)
  }

  input_mat
}

fun_mat_norm <- function(input_mat,method =c("scale","min-max","towards"),axis=1,
                         scale_args = list()){
  stopifnot("input_mat shold be matrix" == is.matrix(input_mat))
  stopifnot("axis only support 0 and 1" == axis %in% c(0,1))
  # 1) choose the methods
  method = match.arg(method,c("scale","min-max","towards"))

  # 2) choose the method
  if(method == "scale"){
    scale_args$x = input_mat
    # run by column
    if(axis == 1){
      message("scale by column")
      input_mat = do.call(scale,scale_args)
    }else{
      # run by row
      message("scale by row")
      input_mat = t(input_mat)
      input_mat = do.call(scale,scale_args)
      input_mat = t(input_mat)
    }

  }else if(method == "min-max"){
    message("normalized into range(0,1),if found record with all same values, the results would be all NaN")
    ft_mm <- function(input_num){
      input_num = input_num - min(input_num,na.rm = T)
      input_num = input_num/max(input_num,na.rm = T)
      input_num
    }
    # run by column
    if(axis == 1){
      message("scale by column")
      for(i in seq_along(input_mat[1,])){
        input_mat[,i] = ft_mm(input_mat[,i])
      }
    }else{
      # run by row
      message("scale by row")
      for(i in seq_along(input_mat[,1])){
        input_mat[i,] = ft_mm(input_mat[i,])
      }
    }

  }else if(method == "towards"){

  }else{
    stop("not found valid input please check input_matrix or method arguments")
  }

  return(input_mat)
}

fun_norm_df <- function(input_df){

}
