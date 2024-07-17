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
fun_mat_norm <- function(input_mat,
                         method =c("scale","min-max","towards"),
                         axis=1,
                         towards_sample = NULL,
                         fun_norm_scale_args = list(),
                         fun_norm_minmax_args = list(),
                         fun_norm_substract_args = list()){
  stopifnot("input_mat shold be matrix" = (is.matrix(input_mat)))
  stopifnot("axis only support 0 and 1" = (axis %in% c(0,1)))
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

#' reanmes genes number into symbols
#'
#' @param input_mat the input matrix
#' @param fun_gene2symbol_args the argument passing tp fun_gene2symbol from biodata
#'
#' @return matrix
#' @export
#'
#' @examples
fun_mat_rename2symbol <- function(input_mat,fun_gene2symbol_args = list()){
  # 1) get the matrix rownames
  matrix_rownames <- rownames(input_mat)
  message("How many genes input: ",length(matrix_rownames))
  fun_gene2symbol_args$input_genes <- matrix_rownames

  # 2) convert the genes into symbol
  suppressMessages({
    gene_map_df <- do.call(biodata::fun_gene2symbol,fun_gene2symbol_args) %>%
      dplyr::filter(!is.na(symbol))
  })
  message("How many records kept: ",nrow(gene_map_df))
  message("How many records kept(unique,final matrix output row length): ",
          length(unique(gene_map_df$symbol)))
  symbols_count <- gene_map_df %>%
    dplyr::count(symbol) %>%
    as.data.frame() %>%
    magrittr::set_rownames(.$symbol)

  # 3) get the symbol uniquely found
  symbol_kept_unique <- gene_map_df %>%
    mutate(symbol_count = symbols_count[.$symbol,"n"]) %>%
    dplyr::filter(symbol_count == 1 | alias_is_symbol)

  mat1 <- input_mat[rownames(symbol_kept_unique),,drop=F]
  rownames(mat1) <- symbol_kept_unique$symbol

  # 4)
  symbol_multiple_parsing <-
    gene_map_df %>%
    dplyr::filter(!symbol %in% symbol_kept_unique$symbol)
  gene_map_info <- apply(input_mat[rownames(symbol_multiple_parsing),],1,sum,na.rm=T)
  message("How many genes with duplicated symbols: ",nrow(symbol_multiple_parsing))
  gene_multi_summary_info <-
    symbol_multiple_parsing %>%
    dplyr::mutate(gene_summary = gene_map_info[rownames(.)]) %>%
    dplyr::group_by(symbol) %>%
    dplyr::arrange(dplyr::desc(gene_summary),.by_group = T) %>%
    dplyr::do(head(.,1)) %>%
    as.data.frame() %>%
    set_rownames(.$alias)
  message("How many symbols are mapped duplicatedly: ",nrow(gene_multi_summary_info))
  mat2 <-  input_mat[rownames(gene_multi_summary_info),,drop=F]
  rownames(mat2) <- gene_multi_summary_info$symbol
  return(rbind(mat1,mat2[,colnames(mat1),drop=F]))
}


#' QN normalization
#'
#' @param input_mat
#' @param qn_args
#'
#' @return matrix
#' @export
#'
#' @examples
fun_mat_qn <- function(input_mat,qn_args = list(keep.names = TRUE)){
  message("perform quantile normalizations")
  qn_args$x = input_mat
  input_mat = do.call(preprocessCore::normalize.quantiles,
                      qn_args)
  return(input_mat)
}
