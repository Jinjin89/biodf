#' Perform the fisher test, return the data.frame
#'
#' @param input_mat the input matrix
#' @param input_name the name of the names(features and rownames of the output)
#'
#' @return data.frame summary of fisher.test results
#' @export
#'

fun_stat_fisher_mat = function(input_mat,input_name){
  # fill na with zero
  input_mat[is.na(input_mat)] = 0
  if (max(dim(input_mat)) > 2){
    hybrid = T
  }else{
    hybrid = F
  }
  res = fisher.test(input_mat,hybrid = hybrid)
  if (max(dim(input_mat)) > 2){
    or = NA
    ci_low = NA
    ci_up = NA
  }else{
    or  = res$estimate
    ci_low = as.numeric(res$conf.int[1])
    ci_up = as.numeric(res$conf.int[2])
  }
  data.frame(row.names = input_name,
             features =input_name,
             or  = or,
             ci_low = ci_low,
             ci_up = ci_up,
             pval = res$p.value)

}

fun_stat_fisher_by_row = function(input_df,
                                  input_abcd = c(2,3,4,5),
                                  input_features = 1){
  # 1) for loop for all the data


}

#' Give a data.frame, compare the variable
#'
#' @param input_df  data.frame
#' @param input_xs discrete variables
#' @param input_ys discrete variables
#' @param echo whether print the confusion matrix
#'
#' @return comparison of the fisher results
#' @export
#'

fun_stat_fisher_by_column = function(input_df,
                                     input_xs,
                                     input_ys,
                                     echo=F){

  # 1) for loop for all the x and y, get pval for each
  purrr::map_df(input_xs,function(x){
    purrr::map_df(input_ys,function(y){
      # 1) get value
      input_mat =
        dplyr::select(input_df,all_of(c(x,y)))
      # 2) remove na value
      input_mat = na.omit(input_mat)
      # 3) get true matix
      input_mat = as.matrix(table(input_mat))

      # 4) get results
      if(echo) print(input_mat)

      f_t = fun_stat_fisher_mat(input_mat,input_name = paste0(x," v.s. ",y))

    })
  })

}
