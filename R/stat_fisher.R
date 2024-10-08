#' Perform the fisher test, return the data.frame
#'
#' @param input_mat the input matrix
#' @param input_name the name of the names(features and rownames of the output)
#' @param chisq whether using chisq.test,default: FALSE
#' @param test_args args passing to fisher.test or chisq.test
#'
#' @return data.frame summary of fisher.test results
#' @export
#'

fun_stat_fisher_mat = function(input_mat,input_name,chisq=F,test_args = list()){
  # fill na with zero
  input_mat[is.na(input_mat)] = 0
  if(chisq){
    message("Using chisq-test")
    #res = stats::chisq.test(input_mat)
    test_args$x = input_mat
    res = do.call(chisq.test,test_args)
    or  = NA
    ci_low = NA
    ci_up = NA
    pval = res$p.value
  }else{
    if (max(dim(input_mat)) > 2) hybrid = T else hybrid = F
    # res = fisher.test(input_mat,hybrid = hybrid)
    test_args$x = input_mat
    test_args$hybrid = hybrid
    res = do.call(fisher.test,test_args)
    if (max(dim(input_mat)) > 2){
      or = NA
      ci_low = NA
      ci_up = NA
      pval =res$p.value
    }else{
      or  = res$estimate
      ci_low = as.numeric(res$conf.int[1])
      ci_up = as.numeric(res$conf.int[2])
      pval =res$p.value
    }
  }

  data.frame(row.names = input_name,
             features =input_name,
             or  = or,
             ci_low = ci_low,
             ci_up = ci_up,
             pval = pval)

}

#' Fisher test by row
#'
#' @param input_df input_data
#' @param input_abcd the abcd column
#' @param chisq whether using chisq.text,defulat: FALSE
#' @param fun_stat_fisher_mat_args fun_stat_fisher_mat_args
#'
#' @return data.frame which appedn the or, ci_low, ci_up, pval column
#' @export
#'
fun_stat_fisher_by_row = function(input_df,
                                  input_abcd = c(2,3,4,5),
                                  chisq=F,
                                  fun_stat_fisher_mat_args = list()){
  input_df$or = NA
  input_df$ci_low = NA
  input_df$ci_up = NA
  input_df$pval = NA

  # 1) for loop for all the data
  for(i in seq_along(input_df[[1]])){
    # 1) get matrix
    a = input_df[[input_abcd[1]]][i]
    b = input_df[[input_abcd[2]]][i]
    c = input_df[[input_abcd[3]]][i]
    d = input_df[[input_abcd[4]]][i]
    mat_tmp = base::matrix(c(a,b,c,d),nrow = 2,ncol = 2,byrow = T)

    # 2) get results
    #fisher_res = fun_stat_fisher_mat(input_mat = mat_tmp,input_name = "tmp")
    fun_stat_fisher_mat_args$input_mat = mat_tmp
    fun_stat_fisher_mat_args$input_name = "tmp"
    fun_stat_fisher_mat_args$chisq = chisq
    fisher_res = do.call(fun_stat_fisher_mat,fun_stat_fisher_mat_args)

    # 3) merge results into input_df
    input_df$or[i] = fisher_res$or
    input_df$ci_low[i] = fisher_res$ci_low
    input_df$ci_up[i] = fisher_res$ci_up
    input_df$pval[i] = fisher_res$pva
  }
  input_df

}

#' Give a data.frame, compare the variable
#'
#' @param input_df  data.frame
#' @param input_xs discrete variables
#' @param input_ys discrete variables
#' @param echo whether print the confusion matrix
#' @param chisq whether using chisq.text,defulat: FALSE
#' @param fun_stat_fisher_mat_args fun_stat_fisher_mat_args
#'
#' @return data.frame,comparison of the fisher results
#' @export
#'

fun_stat_fisher_by_column = function(input_df,
                                     input_xs,
                                     input_ys,
                                     echo=F,
                                     chisq=F,
                                     fun_stat_fisher_mat_args = list()){

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

      #f_t = fun_stat_fisher_mat(input_mat,input_name = paste0(x," v.s. ",y))
      fun_stat_fisher_mat_args$input_mat = input_mat
      fun_stat_fisher_mat_args$input_name = paste0(x," v.s. ",y)
      fun_stat_fisher_mat_args$chisq = chisq
      f_t = do.call(fun_stat_fisher_mat,fun_stat_fisher_mat_args)
      return(f_t)
    })
  })

}
