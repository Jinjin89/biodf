#' Rename the vector the the count of each unique value, output is the same lenght as input
#'
#' @param input_vec the input vector of labels
#' @param count_prefix prefix
#' @param count_suffix suffix
#'
#' @return input vector
#' @export
#'
fun_relabel_vec_count <- \(input_vec,count_prefix = '(',count_suffix = ')'){
  data.frame(data = input_vec) %>%
    dplyr::count(data) %>%
    fun_2df(index_col = 'data') %>%
    dplyr::mutate(newName = paste0(rownames(.),count_prefix,n,count_suffix)) %>%
    dplyr::pull(newName) %>%
    return
}


