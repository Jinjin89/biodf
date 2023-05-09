#' barplot with pvalue
#'
#' @param input_df data.frame
#' @param input_xs x
#' @param input_ys y
#' @param bar_width bar withth
#' @param tn theme
#'
#' @return ggplot obj
#' @export
#'
#' @examples
fun_plot_bar_pval = function(input_df,
                             input_xs,
                             input_ys,
                             palette = "jco",
                             bar_width = 0.3,
                             tn = tn_bar()){
  # 1) barplot for input_xs and input_ys_respectively
  purrr::map(input_xs,function(x){
    purrr::map(input_ys,function(y){
      data_new  = data.frame(x = input_df[[x]],
                             y = input_df[[y]])
      data_new = na.omit(data_new)
      pval = fun_stat_fisher_by_column(data_new,"x","y")
      pval = pval$pval
      x_length = length(unique(data_new$x))
      ggplot2::ggplot(data_new,aes(x,fill = y))+
        ggplot2::geom_bar(width = bar_width,position = "fill")+
        ggplot2::annotate("text",x = (x_length + 1)/2,
                          y = 1.05,
                          label = paste0("p=",format(round(pval,3),nsmall = 3)))+
        tn+
        ggplot2::scale_y_continuous(labels = scales::percent)+
        ggplot2::labs(x = x ,y = y )+
        ggpubr::fill_palette(palette)
    })
  }) -> p_list
  if(length(input_xs) == 1){
    p_list = p_list[[1]]
  }
  p_list
}

#' Volcanoplot
#'
#' @param input_df  input_df
#' @param pval pval select
#' @param pval_cutoff pval cutoff
#' @param logFC logFC select
#' @param logFC_cutoff logFC selected
#' @param palette color
#' @param tn theme default
#'
#' @return ggplot obj
#' @export
#'
#' @examples
fun_plot_volcano1 =function(
    input_df,
    pval = "adj.P.Val",
    pval_cutoff=0.01,
    logFC = "logFC",
    logFC_cutoff= 1,
    palette = c("blue","gray","red"),
    tn = ggplot2::theme_minimal()){

  # cutoff fun
  input_df = fun_deg_cutoff(input_df=input_df,pval = pval,pval_cutoff = pval_cutoff,logFC = logFC,logFC_cutoff = logFC_cutoff)
  # 1) init p
  p = ggplot2::ggplot(input_df,
                      ggplot2::aes(x = !!as.name(logFC),
                          y = -log2(!!as.name(pval))))+
    geom_point(size = 0.7,alpha = 0.2,
               aes(color =p_sig))+
    ggplot2::geom_vline(xintercept = abs(logFC_cutoff),
                        color = "gray",
                        lty = "dashed")+
    ggplot2::geom_vline(xintercept = -abs(logFC_cutoff),
                        color = "gray",
                        lty = "dashed")+
    ggplot2::geom_hline(yintercept = -log2(pval_cutoff),
                        color = "gray",
                        lty = "dashed")+
    tn+
    ggpubr::color_palette(palette)
  p
}



