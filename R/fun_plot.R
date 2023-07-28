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



#' Forestplot
#'
#' @param input_df input_df
#' @param x x
#' @param y y
#' @param ci_low lower confident interval
#' @param ci_up upper confident interval
#' @param y_order the y label order
#' @param errorbar T or F, whether using geom_errorbar function to plot the interval
#' @param segment T or F, whether using geom_segment function to plot the interval
#' @param color_point T or F, whether to color the point
#' @param color_palette the color palete passing to ggpubr::color_palette
#' @param color_use which features used to plot the color
#' @param point_size point size
#' @param point_shape point shape
#' @param vline_color vertical line(x=1)
#'
#' @return ggplot obj
#' @export
#'
#' @examples
fun_plot_forest <-  function(
    input_df,
    x = "hr",
    y = "features",
    ci_low = "ci_low",
    ci_up = "ci_up",
    y_order =NULL,

    errorbar = T,
    segment = F,

    color_point = F,
    color_palette = "default",
    color_use = y,
    point_size= 1,
    point_shape=15,

    vline_color = "gray"


){
  # 1)data preprocess
  if(length(y_order) == 0 ){
    if(!is.factor(input_df[[y]])){
      input_df[[y]] =
        forcats::fct_reorder(
          input_df[[y]],
          input_df[[x]]
        )
    }

  }else{
    y_order = unique(y_order)
    if(length(y_order)!= length(unique(input_df[[y]]))){
      warning("check the y_order params, differ from input_df[[y]]")
    }
    input_df[[y]] =
      factor(input_df[[y]],
             levels = y_order)

  }

  # 2)plot_step
  p = input_df %>%
    ggplot(aes(x =!!as.name(x) ,y =!!as.name(y) ,))+
    geom_vline(xintercept = 1,color = vline_color,lty = "dashed")+
    tn()+
    scale_x_log10()

  # 2.1) whether to plot errobar
  if(errorbar){
    message("Plot using errorbar")
    p = p + geom_errorbar(aes(xmin = !!as.name(ci_low),
                              xmax = !!as.name(ci_up)),
                          width = 0.3)
  }

  # 2.2) whether to plot using segment
  if(segment){
    message("Plot using segment")
    p = p + geom_segment(aes(x = !!as.name(ci_low),
                             xend = !!as.name(ci_up),
                             yend = !!as.name(y)),
                         col  ="gray30")
  }

  # 2.3) point: whether to plot with color
  if(color_point){
    p = p+geom_point(aes(color = !!as.name(y)),size = point_size,shape =point_shape)+
      ggpubr::color_palette(color_palette)+
      theme(legend.position = "none")
  }else{
    p = p + geom_point(size = point_size,shape =point_shape)
  }
  p
}


#' Heatmap plot using Heatmpa in ComplexHeatmap package
#'
#' @param input_matrix the plot matrix
#' @param top_df top_annotaiton data.frame,colnames should be in rownames of input_matrix
#' @param left_df left_annotation data.frame, rownames should be in rownames of input_matrix
#' @param name label the matrix heatmap legend title
#' @param show_column_names whether showing the column name
#' @param show_column_dend whether showing the  column dendro
#' @param cluster_column_slices whether cluster column slice
#' @param show_row_dend show_row_dend
#' @param row_title row_title, defualt is " ", showing nothing
#' @param cluster_row_slices cluster_row_slices
#' @param row_name_cex rownames font params, parssing to gp()
#' @param row_name_face rownames font params, parssing to gp()
#' @param row_name_pval_cutoff pvalue cutoff if the row_name was calculated
#' @param row_name_color_sig_up red
#' @param row_name_color_ns black
#' @param row_name_color_sig_down blue
#' @param top_split  left split data, should be in colnames of left_df or number index
#' @param top_col list of top_annotation color map
#' @param top_args params list passing to HeatmapAnnotation
#' @param left_split left split data, should be in colnames of left_df or number index
#' @param left_col list of left_annotataion color map
#' @param left_args params list passing to RowAnnotation
#' @param ... passing to ComplexHeatmap::Heatmap
#' @param hide_column_df_annotation F
#' @param hide_row_df_annotation F
#'
#' @return ComplexHeatmap plot
#' @export
#'
#' @examples
fun_plot_heatmap = function(
    input_matrix,
    top_df,
    left_df,

    name = "z-score",
    show_column_names=F,
    show_column_dend=F,
    cluster_column_slices =F,

    show_row_dend=F,
    row_title=" ",
    cluster_row_slices =F,

    row_name_cex = 0.5,
    row_name_face="italic",
    row_name_pval_cutoff= 0.05,
    row_name_color_sig_up = "red",
    row_name_color_ns = "black",
    row_name_color_sig_down = "blue",

    top_split=1,
    top_col = biodf::fun_utils_get_color(top_df),
    top_args = list(),


    left_split=1,
    left_col = biodf::fun_utils_get_color(left_df),
    left_args = list(),

    hide_column_df_annotation = F,
    hide_row_df_annotation = F,
    ...

){
  suppressMessages(require(ComplexHeatmap))
  # 1) get matrix
  top_names = intersect(colnames(input_matrix),
                        rownames(top_df))
  left_names = intersect(rownames(input_matrix),
                         rownames(left_df))
  input_matrix = input_matrix[left_names,top_names]
  top_df = top_df[top_names,,drop =F]
  left_df = left_df[left_names,,drop =F]


  # 2) get args
  # 2.1) top args
  top_args$df = top_df
  top_args$col = top_col
  row_name_color =row_name_color_ns
  if(length(top_split) != 0 ){
    top_split = top_df[[top_split]]
    top_df_for_compare_means =
      data.frame(sample = rownames(top_df),
                 feature = top_split) %>%
      biodf::fun_merge_matrix(
        as.matrix(input_matrix),
        input_variables = rownames(left_df)
      ) %>%
      biodf::fun_stat_compare_means(
        input_variables = rownames(left_df),
        input_groups = "feature") %>%
      dplyr::mutate(p_sig = ifelse(pval<row_name_pval_cutoff,"*","N.S."))

    # parsing color, if the feature is of two factor, get up and down, others, no
    group_labels = as.character(sort(unique(top_split)))
    if(length(group_labels) == 2){
      # 1) get name
      top_df_for_compare_means$direction = top_df_for_compare_means[[group_labels[2]]] -top_df_for_compare_means[[group_labels[1]]]
      top_df_for_compare_means$direction = ifelse(top_df_for_compare_means$direction > 0,"up","down")

      top_df_for_compare_means %<>%
        dplyr::mutate(color = dplyr::case_when(
          p_sig == "N.S."~row_name_color_ns,
          direction == "up"~row_name_color_sig_up,
          direction == "down"~row_name_color_sig_down,
          TRUE~NA_character_
        ))
      row_name_color = top_df_for_compare_means[left_names,"color"]
    }else{
      message("not calculate the up or down")
      top_df_for_compare_means$color = ifelse(top_df_for_compare_means$p_sig == "*",
                                              row_name_color_sig_up,
                                              row_name_color_ns)
      row_name_color = top_df_for_compare_means[left_names,"color"]
    }

  }

  # 2.2) left args
  left_args$df = left_df
  left_args$col = left_col
  if(length(left_split) != 0 ){
    left_split = left_df[[left_split]]

  }

  # 2.3) if hide annotations

  if(hide_column_df_annotation){
    top_args$df = NULL
  }
  if(hide_row_df_annotation){
    left_args$df = NULL
  }


  # 3)plot_step
  ComplexHeatmap::Heatmap(

    # 1) main
    matrix = input_matrix,
    name = name,
    show_column_names = show_column_names,
    show_column_dend = show_column_dend,
    show_row_dend = show_row_dend,


    # 2) top
    top_annotation = tryCatch(do.call("HeatmapAnnotation",top_args),error = function(e) NULL),
    column_split =top_split,
    cluster_column_slices =cluster_column_slices,


    # 3) left_annotation
    left_annotation = tryCatch(do.call("rowAnnotation",left_args),error = function(e) NULL),
    row_split = left_split,
    row_title =row_title,
    cluster_row_slices = cluster_row_slices,
    row_names_gp = gpar(cex = row_name_cex,
                        fontface = row_name_face,
                        col = row_name_color),
    ...
  )

}


#' boxplot
#'
#' @param input_df data.frame
#' @param input_x group data
#' @param input_y ydata,regular expression
#' @param palette palette passing to fill_paltte in ggpubr
#' @param theme_now ggplot theme
#' @param str_remove ggtitle remove
#' @param p_text_size text size
#'
#' @return ggplot or list with ggplots
#' @export
#'
#' @examples
fun_plot_boxplot <-function(input_df,input_x,input_y,palette="jco",
                            theme_now = tn(),
                            remove_suffix=T,
                            str_remove = "\\$|_",p_text_size = 3){
  # 1) loop x
  purrr::map(input_x,function(each_x){
    # 2) loop y
    purrr::map(input_y,function(each_y){
      plot_tmp =
        input_df %>%
        dplyr::select(all_of(each_x),matches(each_y)) %>%
        tidyr::pivot_longer(-!!as.name(each_x))

      if(remove_suffix){
        plot_tmp$name = stringr::str_remove(plot_tmp$name,each_y)
      }
      plot_tmp %>%
        ggplot(aes(name,value,fill =!!as.name(each_x)))+
        geom_boxplot(outlier.size = 0.5)+
        theme_now+
        theme(axis.text.x = element_text(angle = 75,hjust = 1))+
        ggpubr::fill_palette(palette)+
        ggtitle(stringr::str_remove_all(each_y,str_remove))+
        ggpubr::stat_compare_means(label = "p.signif",size = p_text_size)
    }) %>%
      purrr::set_names(input_y)

  }) %>%
    purrr::set_names(input_x) -> p_list

  if(length(input_y) == 1){
    for(i in seq_along(input_x)){
      p_list[[i]] = p_list[[i]][[1]]
    }

  }

  if(length(input_x) == 1){
    p_list <- p_list[[1]]
  }
  p_list
}
