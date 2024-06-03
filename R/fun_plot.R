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
                             tn = tn_bar(),
                             chisq=F,
                             fun_stat_fisher_by_column_args = list()){
  # 1) barplot for input_xs and input_ys_respectively
  purrr::map(input_xs,function(x){
    purrr::map(input_ys,function(y){
      data_new  = data.frame(x = input_df[[x]],
                             y = input_df[[y]])
      data_new = na.omit(data_new)
      # pval = fun_stat_fisher_by_column(data_new,"x","y")
      fun_stat_fisher_by_column_args$input_df = data_new
      fun_stat_fisher_by_column_args$input_xs = "x"
      fun_stat_fisher_by_column_args$input_ys = "y"
      fun_stat_fisher_by_column_args$chisq = chisq
      pval = do.call(fun_stat_fisher_by_column,
                     fun_stat_fisher_by_column_args)
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
  if(length(input_ys) == 1){
    p_list = unlist(p_list,recursive = F)
  }
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





#' Heatmap plot with fill and p-value text
#'
#' @param input_cor_df correlation reults
#' @param x x
#' @param y y
#' @param fill correlation
#' @param text pval
#' @param Heatmap_args other params passing to args
#' @param Heatmap_args_additional additional args into heatmap
#'
#' @return heatmap plot
#' @export
#'
#' @examples
fun_plot_heatmap_tri <- function(
    input_cor_df,x = "x",
    y = "y",fill = "cor",text="p_sig",
    Heatmap_args = list(
      cluster_rows = F,
      cluster_columns = F,
      name = "correlation",
      row_names_gp = gpar(fontsize = 8,fontface = "italic"),
      column_names_gp = gpar(fontsize = 8,fontface = "plain"),
      rect_gp = gpar(fill = "white",color = "gray"),
      col = circlize::colorRamp2(
        c(-0.5,0,0.5),colors = c(scales::muted("blue"),"white",scales::muted("red")))),
    Heatmap_args_additional = list(),
    hm_text_size = 7){

  suppressMessages(require(circlize))
  suppressMessages(require(ComplexHeatmap))

  cor_matrix <-
    input_cor_df %>%
    pivot_wider(id_cols = !!as.name(y),
                names_from = !!as.name(x),
                values_from = !!as.name(fill)) %>%
    as.data.frame() %>%
    magrittr::set_rownames(.[[y]]) %>%
    dplyr::select(-!!as.name(y))

  # 3) p value text matrix
  p_matrix <-
    input_cor_df %>%
    pivot_wider(id_cols = !!as.name(y),
                names_from = !!as.name(x),
                values_from = !!as.name(text)) %>%
    as.data.frame() %>%
    magrittr::set_rownames(.[[y]]) %>%
    dplyr::select(-!!as.name(y))

  # 4) rownames
  col_labels = colnames(cor_matrix) %>% sort
  row_labels = rownames(cor_matrix) %>% sort
  cor_matrix = cor_matrix[row_labels,col_labels]
  p_matrix = p_matrix[row_labels,col_labels]

  if(length(names(Heatmap_args_additional)) > 0){
    for(each_arg in names(Heatmap_args_additional)){
      message("Adding args: ", each_arg)
      Heatmap_args[[each_arg]] = Heatmap_args_additional[[each_arg]]
    }
  }
  Heatmap_args$matrix = as.matrix(cor_matrix)


  if(is.null(Heatmap_args$cell_fun)){
    Heatmap_args$cell_fun <- function(j, i, x, y, w, h, fill) {
      grid.polygon(
        x = unit.c(x - 0.5*w, x - 0.5*w, x + 0.5*w),
        y = unit.c(y - 0.5*h, y + 0.5*h, y + 0.5*h),
        gp = gpar(
          col = "gray",
          fill = Heatmap_args$col(cor_matrix[i,j])
        ))

      grid.text(
        label = p_matrix[i,j],
        x = unit.c(x + 0.5*w),
        y = unit.c(y - 0.5*h),
        gp = gpar(fontsize = hm_text_size),
        hjust = 1.2,vjust = 0)

      grid.text(
        label = round(cor_matrix[i,j],2),
        x = unit.c(x - 0.5*w),
        y = unit.c(y + 0.5*h),
        gp = gpar(fontsize = hm_text_size),
        hjust = 0,
        vjust = 1.2)
    }
  }
  do.call(ComplexHeatmap::Heatmap,
          Heatmap_args)

}

fun_get_legend <- function(plot){
  if (inherits(plot, 'gg')) {
    gt <- ggplot_gtable(ggplot_build(plot))
  } else {
    gt <- ggplotify::as.grob(plot)
  }
  gname <- vapply(gt$grobs, function(x) x$name, FUN.VALUE = character(1))
  idx <- which(gname == "guide-box")
  legend <- gt$grobs[[idx]]
  return(legend)
}

#' Butterfly plot for correlation anaysis
#'
#' @param input_df data.frame or correlation results from `fun_stat_cor`
#' @param input_key key
#' @param input_variables_lower lower heatmap variables
#' @param input_variables_upper upper heatmap variables
#' @param input_df_is_cor input_df_is_cor
#' @param return_plot_list whethre return plot list
#' @param plot_text_size plot_text_size
#' @param plot_x_angle plot_x_angle
#' @param plot_color_map_bre plot_color_map_bre
#' @param plot_color_map_col plot_color_map_col
#' @param tn_lgd tn_lgd
#' @param plot_rect_col plot_rect_col
#' @param plot_curvature plot_curvature
#' @param lower_pos_shift lower_pos_shift
#' @param upper_pos_shift upper_pos_shift
#' @param larger_length larger_length
#' @param smaller_length smaller_length
#'
#' @return plot
#' @export
#'
#' @examples

fun_plot_butterfly <- \(input_df,input_key,input_variables_lower,input_variables_upper,
                        # plot params
                        return_plot_list=F,
                        plot_text_size = 7,
                        plot_x_angle = 75,
                        plot_color_map_bre = c(-1,0,1),
                        plot_color_map_col = c("#3A3A98","white","#832424"),
                        tn_lgd = theme(legend.key.size = unit(5,"mm"),
                                       legend.position = "bottom",
                                       legend.text = element_text(size = 9),
                                       legend.title = element_text(size = 9),
                                       legend.spacing = unit(5,"mm"),
                                       legend.direction = "vertical"),
                        # rect_color
                        plot_rect_col = "gray",
                        plot_curvature = 0.1,
                        tn = tn_empty(),
                        # arrange
                        lower_pos_shift = -2,
                        upper_pos_shift = 2,
                        larger_length = 3,
                        smaller_length = 1.6,

                        input_df_is_cor=F){
  # loading required packages
  require(grid)
  input_variables_lower = unique(input_variables_lower)
  input_variables_upper = unique(input_variables_upper)
  all_var <- c(input_key,input_variables_lower,input_variables_upper)
  lower_var_count = length(input_variables_lower)
  upper_var_count = length(input_variables_upper)
  total_var_count = length(all_var)

  #########correlation data parsing########
  if(input_df_is_cor){
    message("The input is results is correlation,")
    cor_res <- input_df
  }else{
    message("The input is data.frame, perform correlation")
    suppressWarnings(
      {
        cor_res <- input_df %>% fun_stat_cor(all_var,all_var)
      }
    )

  }
  plot_color_map <- circlize::colorRamp2(plot_color_map_bre,plot_color_map_col)

  # 2) subset the data into lower part
  #########subsetting data########
  cor_res_lower <- cor_res %>%
    dplyr::filter(x %in% c(input_key,input_variables_lower)) %>%
    dplyr::filter(y %in% c(input_key,input_variables_lower)) %>%
    mutate(fill = plot_color_map(.$cor))

  cor_res_upper <- cor_res %>%
    dplyr::filter(x %in% c(input_key,input_variables_upper)) %>%
    dplyr::filter(y %in% c(input_key,input_variables_upper))%>%
    mutate(fill = plot_color_map(.$cor))

  message("Change input var into factors")
  cor_res_lower$x = factor(cor_res_lower$x,levels = c(input_variables_lower,input_key))
  cor_res_lower$y = factor(cor_res_lower$y,levels = c(input_key,input_variables_lower[lower_var_count:1]))

  cor_res_upper$x = factor(cor_res_upper$x,levels = c(input_variables_upper,input_key))
  cor_res_upper$y = factor(cor_res_upper$y,levels = c(input_key,input_variables_upper[upper_var_count:1]))

  #########lower data sorting########
  message("Sorting lower-plot data")
  cor_res_lower %<>%
    dplyr::mutate(x_number = as.numeric(x)) %>%
    dplyr::mutate(y_number = as.numeric(y)) %>%
    dplyr::filter(y_number+x_number <= lower_var_count+2)

  lower_point <-
    cor_res_lower %>%
    dplyr::mutate(point_x = (x_number+1),
                  point_y = lower_var_count-x_number+2) %>%
    dplyr::filter(str_detect(rownames(.),paste0(input_key))) %>%
    dplyr::select(point_x,point_y)%>%
    rbind(data.frame(
      point_x = lower_var_count + lower_pos_shift+1,
      point_y = lower_var_count + lower_pos_shift+1
    ))

  lower_line = cor_res_lower %>%
    dplyr::mutate(point_x = (x_number+1),
                  point_y = lower_var_count-x_number+2) %>%
    dplyr::filter(str_detect(rownames(.),paste0(input_key))) %>%
    dplyr::mutate(key_x = lower_var_count + lower_pos_shift+1,
                  key_y = lower_var_count + lower_pos_shift+1) %>%
    dplyr::mutate(relation = ifelse(.$cor >0 ,"Positive","Negtive"),
                  correlation = abs(.$cor),
                  pvalue = fun_utils_p2star(.$pval,"N.S.")) %>%
    dplyr::mutate(relation = factor(.$relation,levels = c("Positive","Negtive"))) %>%
    dplyr::mutate(pvalue = factor(pvalue,levels = c("N.S.","*","**","***","****")))

  #########upper data sorting########
  message("Sorting upper-plot data")
  cor_res_upper %<>%
    dplyr::mutate(x_number = as.numeric(x)) %>%
    dplyr::mutate(y_number = as.numeric(y)) %>%
    dplyr::filter(y_number + x_number >= upper_var_count+2)

  upper_point <-
    cor_res_upper %>%
    dplyr::mutate(point_x = (y_number -1),
                  point_y = upper_var_count - y_number+2) %>%
    dplyr::filter(str_detect(rownames(.),paste0(input_key))) %>%
    dplyr::select(point_x,point_y) %>%
    rbind(data.frame(
      point_x = upper_pos_shift,
      point_y = upper_pos_shift
    ))
  upper_line = cor_res_upper %>%
    dplyr::mutate(point_x = upper_var_count - y_number+1,
                  point_y = y_number) %>%
    dplyr::filter(str_detect(rownames(.),paste0(input_key))) %>%
    dplyr::mutate(key_x = upper_pos_shift,
                  key_y = upper_pos_shift) %>%
    dplyr::mutate(relation = ifelse(.$cor >0 ,"Positive","Negtive"),
                  correlation = abs(.$cor),
                  pvalue = fun_utils_p2star(.$pval,"N.S.")) %>%
    dplyr::mutate(relation = factor(.$relation,levels = c("Positive","Negtive"))) %>%
    dplyr::mutate(pvalue = factor(pvalue,levels = c("N.S.","*","**","***","****")))
  #########plot lower########

  p_lower <-
    cor_res_lower %>%
    ggplot(aes(x,y))+
    geom_tile(aes(fill = fill),color=plot_rect_col)+
    scale_fill_identity()+
    scale_x_discrete(expand = c(0,0))+
    scale_y_discrete(expand = c(0,0))+
    geom_curve(
      data = lower_line,
      aes(x = key_x,y=key_y,xend = point_x,yend = point_y,
          lty = relation,linewidth = correlation,color = pvalue),
      curvature = plot_curvature,
      inherit.aes = F,show.legend = F
    )+
    scale_linewidth(range = c(1,2))+
    geom_point(data = lower_point,aes(point_x,point_y),
               inherit.aes = F)+
    tn+
    theme(axis.text.x = element_text(angle = plot_x_angle,size = plot_text_size,hjust = 1),
          axis.text.y = element_text(size = plot_text_size,hjust = 1))+
    coord_fixed(clip = "off")
  #########plot upper########

  p_upper <-
    cor_res_upper %>%
    ggplot(aes(x,y))+
    geom_tile(aes(fill = fill),color=plot_rect_col)+
    scale_x_discrete(position = "top",expand = c(0,0))+
    scale_y_discrete(position = "right",expand = c(0,0))+
    scale_fill_identity()+
    coord_fixed(clip = "off")+
    geom_curve(
      data = upper_line,
      aes(x = key_x,y=key_y,xend = point_x,yend = point_y,
          lty = relation,linewidth = correlation,color = pvalue),
      curvature = plot_curvature,
      inherit.aes = F,show.legend = F
    )+
    scale_linewidth(range = c(1,2))+
    geom_point(data = upper_point,aes(point_x,point_y),
               inherit.aes = F)+
    tn+
    theme(axis.text.x = element_text(angle = plot_x_angle,size = plot_text_size,hjust = 0),
          axis.text.y = element_text(size = plot_text_size,hjust = 1))

  #########legend plot########
  legend_df <- data.frame(
    cor = seq(min(plot_color_map_bre),max(plot_color_map_bre),length.out = 20),
    pval = c(0.0001,0.001,0.01,0.05,0.6,rep(1,15))
  ) %>%
    dplyr::mutate(relation = ifelse(.$cor >0 ,"Positive","Negtive"),
                  correlation = abs(.$cor),
                  pvalue = fun_utils_p2star(.$pval,"N.S.")) %>%
    dplyr::mutate(relation = factor(.$relation,levels = c("Positive","Negtive"))) %>%
    dplyr::mutate(pvalue = factor(pvalue,levels = c("N.S.","*","**","***","****"))) %>%
    dplyr::mutate(x =1:nrow(.), y =1:nrow(.))

  p_hm <- legend_df %>%
    ggplot(aes(x,y,xend = x +1, yend = y + 1))+
    geom_tile(aes(fill = cor))+
    scale_fill_gradient2(low = plot_color_map_col[1],
                         mid = plot_color_map_col[2],
                         high = plot_color_map_col[3])+
    theme_void()+
    tn_lgd

  p_line <-
    legend_df %>%
    ggplot(aes(x,y,xend = x +1, yend = y + 1))+
    geom_curve(aes(size = correlation,color = pvalue,lty = relation))+
    theme_void()+
    tn_lgd+    scale_size(range = c(1,2))

  p_hm_lgd = fun_get_legend(p_hm)
  p_line_lgd = fun_get_legend(p_line)

  p_list <- list(l = p_lower,u = p_upper,lgd = list(p_hm_lgd,p_line_lgd))
  #########return##########
  if(return_plot_list){
    p_list
  }else{
    vpp <- function(x,y) viewport(layout.pos.row = x,layout.pos.col = y)
    grid.newpage()
    pushViewport(
      viewport(layout = grid.layout(
        nrow = 2,ncol = 2,
        widths = unit(c(larger_length,smaller_length),"null"),
        heights = unit(c(smaller_length,larger_length),"null"))))
    print(p_list[[1]],vp = vpp(2,1))

    upViewport()
    pushViewport(
      viewport(layout = grid.layout(
        nrow = 2,ncol = 2,
        widths = unit(c(smaller_length,larger_length),"null"),
        heights = unit(c(larger_length,smaller_length),"null"))))
    print(p_list[[2]],vp = vpp(1,2))


    pushViewport(
      viewport(layout = grid.layout(
        nrow = 2,ncol = 2,
        widths = unit(c(smaller_length,larger_length),"null"),
        heights = unit(c(smaller_length,larger_length),"null"))))

    pushViewport(viewport(layout.pos.col = 1,layout.pos.row = 1))
    grid.draw(p_line_lgd)
    upViewport()
    pushViewport(viewport(layout.pos.col = 2,layout.pos.row = 2))
    grid.draw(p_hm_lgd)
    upViewport()
  }
}



#' boxplot and density
#'
#' @param input_df data.frame
#' @param x x
#' @param y y
#' @param palette "jco"
#' @param heights ggarange heights
#' @param tn_top top plot theme
#' @param tn_bottom bottom plot theme
#'
#' @return plot
#' @export
#'
#' @examples
fun_plot_density_box <- \(
  input_df,
  x = "Group",
  y = "Risk_score",
  palette = "jco",
  heights = c(2,1),
  tn_top = tn(),
  tn_bottom = tn()
){

  input_df %>%
    ggplot2::ggplot(
      ggplot2::aes(
        !!as.name(y),
        group = !!as.name(x),
        fill =  !!as.name(x)))+
    ggplot2::geom_density(alpha = 0.2)+
    tn_top+
    ggplot2::theme(legend.position = "top")+
    ggpubr::fill_palette(palette) -> p_d

  input_df %>%
    ggplot2::ggplot(
      ggplot2::aes(
        !!as.name(x),
        y = !!as.name(y),
        fill =  !!as.name(x)))+
    ggplot2::geom_boxplot(width = 0.4)+
    tn_bottom+
    tn_no_legend()+
    ggplot2::coord_flip()+
    ggpubr::fill_palette(palette)+
    ggpubr::stat_compare_means(comparisons = list(c(1,2)),label = "p.signif") -> p_box
  ggpubr::ggarrange(
    p_d,
    p_box,
    nrow = 2,ncol = 1,heights = heights)
}




#' forest plot
#'
#' @param input_df input_df
#' @param x odds ratio or  hr
#' @param y genes
#' @param ci_low ci.low
#' @param ci_up ci.up
#' @param feature1 features
#' @param feature2 features
#' @param pval pval
#' @param order order
#' @param desc T or F
#' @param tn theme
#' @param widths the widths of ggarrange
#'
#' @return ggarrange obj
#' @export
#'

fun_plot_forest2 <- function(
    input_df,
    x = "or",
    y = 'Hugo_Symbol',
    ci_low = "ci.low",
    ci_up = "ci.up",
    feature1 = 'High',
    feature2 = "Low",
    pval = "pval",
    order = NULL,
    desc=F,
    tn = tn_empty(),
    widths = c(2,1,3,1.5)
){

  if(len(order) == 0){
    input_df[[y]] <- forcats::fct_reorder(input_df[[y]],input_df[[x]],.desc = desc)
  }else{
    input_df[[y]] <- factor(input_df[[y]],levels = order)
  }
  tn_no_y <-  theme(axis.text.y = element_blank())
  tn_no_x <-  theme(axis.text.x = element_blank())
  tn_title_mid =theme(plot.title = element_text(hjust = 0.5))

  input_df[[pval]] = paste0("p=",format(input_df[[pval]],digits  = 1,trim = T))

  p_label1 <- input_df %>% ggplot(aes(x=1,y = !!as.name(y)))+tn_empty()+tn_no_x+tn_title_mid
  p_label1 <- p_label1+geom_text(aes(label = !!as.name(feature1)))+
    ggtitle(feature1)

  p_label2 <- input_df %>% ggplot(aes(x=1,y = !!as.name(y)))+tn_empty()+tn_no_y+tn_no_x+tn_title_mid
  p_label2 <- p_label2+geom_text(aes(label = !!as.name(feature2)))+
    ggtitle(feature2)

  p_forest <- input_df %>% ggplot(aes(y = !!as.name(y)))+tn_empty()+tn_no_y+
    geom_vline(xintercept = 1,lty = 'dashed',color = "gray")+tn_title_mid

  p_forest <- p_forest + geom_point(aes(x = !!as.name(
    x)))
  p_forest <- p_forest +
    geom_errorbar(aes(xmin = !!as.name(ci_low),xmax = !!as.name(ci_up)),
                  width = 0.3)+
    scale_x_log10()+
    ggtitle(x)

  p_pval <- input_df %>% ggplot(aes(x=1,y = !!as.name(y)))+tn_empty()+tn_no_y+tn_no_x
  p_pval <- p_pval+geom_text(aes(label = !!as.name(pval)))+
    ggtitle(pval)+tn_title_mid
  ggarrange(p_label1,p_label2,p_forest,p_pval,nrow = 1,ncol = 4,
            widths = widths,align = 'h')

}

