fun_norm_matrix <- function(input_mat,method =c("scale","min-max","towards"),axis=1,
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

