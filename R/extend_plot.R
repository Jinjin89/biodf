#' Plot heatmap with pvalue along the side
#'
#' @param input_df data.frame with x and rownames of input-y
#' @param input_x x
#' @param input_y data.frame, rownames should be in input_df
#' @param method sp or p
#' @param fontsize fontsize
#' @param row_names_gp row_names_gp
#' @param row_split which data in input_y to split row
#' @param left_col left annotation color
#' @param hm_col hm_col map
#'
#' @return NULL
#' @export
#'
#' @examples
fun_extend_plot_heatmap_cor <- function(
    input_df,
    input_x,
    input_y,
    method=  "sp",
    fontsize = 5,
    row_names_gp = gpar(cex = 0.5,fontface = "italic"),
    row_split = 1,
    left_col = NULL,
    hm_col = circlize::colorRamp2(breaks = c(-1,0,1),colors = c("blue","white","red"))){


  # 1) get correlation
  cor_res = fun_stat_cor(input_df,input_x,rownames(input_y),method = method) %>%
    mutate(p_sig = fun_utils_p2star(pval))
  rownames(cor_res) = cor_res$y
  print(cor_res)


  # 2) hm
  # 2.1) hm correlation
  hm_cor = cor_res %>% dplyr::select(cor)
  colnames(hm_cor) = "correlation"
  hm_cor %<>%
    dplyr::arrange(desc(correlation))

  # 2.2) hm pvalues
  hm_pval = cor_res %>% dplyr::select(p_sig)
  colnames(hm_pval) = "pvalue"

  # 2.3) reorder
  hm_pval = hm_pval[rownames(hm_cor),,drop=F]
  input_y  = input_y[rownames(hm_cor),,drop=F]

  # 2.4) get row params
  if(ncol(input_y) == 0){
    left_anno = NULL
  }else{
    if(length(left_col) == 0){
      left_col = fun_utils_get_color(input_y)
    }
    left_anno = rowAnnotation(df = input_y,col =left_col)
  }

  if(length(row_split) == 0){
    row_split = split
  }else{
    row_split = input_y[[row_split]]
  }

  p1 = Heatmap(
    matrix = hm_cor,
    name = "correlation",
    column_title = " ",
    col = hm_col,
    show_column_names = F,
    cluster_rows = F,
    rect_gp = gpar(col = 'gray'),
    cell_fun = function(j, i, x, y, w, h, col) {
      grid.text(round(hm_cor[i, j],2), x, y,gp = gpar(fontsize = fontsize))
    },

    # row annotaiton
    left_annotation = left_anno,
    row_split = row_split
  )

  p2 = Heatmap(
    matrix = hm_pval,
    column_title = " ",
    show_column_names = F,
    cluster_rows = F,
    rect_gp = gpar(col = 'gray',fill = NA),
    cell_fun = function(j, i, x, y, w, h, col) {
      grid.text(hm_pval[i, j], x, y,gp = gpar(fontsize = fontsize))
    },
    show_heatmap_legend = F,
    row_names_gp = row_names_gp,
    row_split = row_split
  )

  draw(p1+p2)

}
