
#' Perform statistical test for variables and each group
#'
#' @param input_df input data.frame
#' @param input_variables input variables
#' @param input_groups input groups
#' @param fun test function, like `t.test`, `wilcox.test`,`kruskal.test`
#' @param fun_sum how to summarise data for each group
#'
#' @return data.frame with pvalue, and summarise each group parms
#' @export
#'
#' @examples
fun_stat_compare_means = function(input_df,input_variables,input_groups="Group",fun = "kruskal.test",fun_sum = "median"){
  # 1) for loop for all the groups
  purrr::map(input_groups,function(each_group){
    # 2) then loop the variables, calculate the variable in each group
    purrr::map(input_variables,function(each_val){
      # 1) get new data
      data_tmp =
        na.omit(data.frame(
          variable = input_df[[each_val]],
          group = input_df[[each_group]]
        ))

      # 2) call with
      test_res = do.call(
        fun,
        args = list(
          formula = variable~group,
          data = data_tmp
        ))

      # 3) summarse the results
      res_df = data.frame(row.names = each_val,features = each_val)
      res_df$pval = test_res$p.value
      group_total = sort(unique(data_tmp$group))
      for(g in group_total){
        res_df[[paste0(g)]] = do.call(
          fun_sum,
          args = list(data_tmp$variable[data_tmp$group == g]))

      }
      res_df$summarise_methods = fun_sum
      res_df$test_methods = fun
      res_df
    }) %>%
      bind_rows()

  }) -> test_results
  names(test_results) = input_groups

  if(length(input_groups) == 1){
    test_results = test_results[[1]]
  }

  test_results
}




#' Calculate the correlation between the input x and y one by one
#'
#' @param input_df the input data.frame
#' @param input_xs input x ,should be in input_df
#' @param input_ys input y ,should be in input_df
#' @param method spearman or person
#'
#' @return data.frame of the correlation results
#' @export
#'
#' @examples
fun_stat_cor = function(input_df,
                        input_xs,
                        input_ys,
                        method = "sp"){
  # 1) get x
  x = dplyr::select(input_df,dplyr::one_of(input_xs))
  y = dplyr::select(input_df,dplyr::one_of(input_ys))

  # 2) get y
  purrr::map_df(colnames(x),function(each_x){
    purrr::map_df(colnames(y),function(each_y){
      res =
        cor.test(
          input_df[[each_x]],
          input_df[[each_y]],
          use = "pairwise.complete.obs",
          method = method)

      data.frame(
        row.names = paste0(each_x,"&",each_y),
        x = each_x,
        y = each_y,
        cor  = res$estimate,
        pval = res$p.value)
    })
  })
}


