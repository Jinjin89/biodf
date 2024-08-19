
#' fisher test
#'
#' @param input_df input_df
#' @param input_variables input_variables
#' @param input_group input_group
#' @param output_variables_names outputnames
#' @param input_pct input_pct
#' @param fun_stat_fisher_by_column_args fun_stat_fisher_by_column_args
#' @param return_grouped_df return_grouped_df
#' @param fun_surv_cutoff_args fun_surv_cutoff_args
#'
#' @return data.frame,results or grouped df
#' @export
#'


fun_df_cutoff_fisher <- function (input_df, input_variables,
                                  input_group = "Group",
                                  output_variables_names =input_variables[1],
                                  return_grouped_df = T,
                                  fun_surv_cutoff_args = list(),
                                  input_pct = seq(0.25,0.75, 0.01),
                                  fun_stat_fisher_by_column_args = list())
{
  df_raw = input_df
  input_variables = input_variables[1]
  input_df <- input_df %>% dplyr::select(dplyr::all_of(c(input_variables,
                                                         input_group))) %>% na.omit()
  fun_df_cutoff_fisher_each <- function(each_pct, df_used = input_df,
                                        stat_args = fun_stat_fisher_by_column_args, group_used = input_group,
                                        var_used = input_variables) {
    data_tmp = df_used[[var_used]]
    data_cutoff = as.numeric(quantile(data_tmp, each_pct,
                                      na.rm = T))
    df_used[[var_used]] <- ifelse(data_tmp < data_cutoff,
                                  "Low", "High")
    df_used[[var_used]] <- factor(df_used[[var_used]], levels = c("Low",
                                                                  "High"))
    fun_stat_fisher_by_column_args$input_df = df_used
    fun_stat_fisher_by_column_args$input_xs = group_used
    fun_stat_fisher_by_column_args$input_ys = var_used
    do.call(fun_stat_fisher_by_column, stat_args) %>% dplyr::mutate(pct = each_pct)
  }
  res <- furrr::future_map(input_pct, fun_df_cutoff_fisher_each) %>%
    purrr::list_rbind() %>% tibble() %>% as.data.frame()
  if (return_grouped_df) {
    pct_select <- res$pct[res$pval == min(res$pval, na.rm = T)]
    if (len(pct_select) > 1) {
      message("multiple pct found(use the first): ",
              paste0(pct_select,collapse = ","))
    }
    message("Optimal cutoff: ",
            pct_select = pct_select[1],
            ", pval is: ", min(res$pval))
    fun_surv_cutoff_args$input_df = df_raw
    fun_surv_cutoff_args$input_variables = input_variables
    fun_surv_cutoff_args$input_pct = pct_select[1]
    fun_surv_cutoff_args$output_variables_names = output_variables_names
    do.call(fun_surv_cutoff, fun_surv_cutoff_args)
  }
  else {
    res
  }
}
