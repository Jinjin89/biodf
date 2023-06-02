#' This is broadcast function,change the shape of input_targets so that could fit into input_variables
#'
#' @param input_variables The input variables
#' @param input_targets The output with same shape as input variables
#'
#' @return return a vector of the input_variables with same shape
#' @export
#'
fun_utils_broadcast = function(input_variables,input_targets){

  if(is.null(input_targets)){
    return(input_variables)
  }else if(length(input_variables) == length(input_targets)){
    return(input_targets)
  }else if(length(input_targets) == 1){
    return(rep(input_variables,length(input_variables)))
  }

}

#' Default theme
#'
#' @param ... passing to theme
#'
#' @return theme object
#' @export


tn = function(...){
  ggplot2::theme_minimal()+
    ggplot2::theme(
      panel.grid = ggplot2::element_blank(),
      axis.line = ggplot2::element_line(),
      axis.ticks = ggplot2::element_line())+
    ggplot2::theme(...)
}

#' theme for barplot
#'
#' @param ... passing to theme
#' @param show_axis_title  whether plot axis title
#' @param legend_text  legend text size
#'
#' @return theme obj
#' @export
#'
tn_bar = function(show_axis_title=F,legend_text = 6,...){
  ggplot2::theme_minimal()+
    ggplot2::theme(
      panel.grid = ggplot2::element_blank(),
      panel.background = ggplot2::element_rect(fill =NA,color = "black"),
      axis.ticks = ggplot2::element_line(color = "black"),
      axis.title = element_blank(),
      legend.title = ggplot2::element_blank(),
      legend.text = ggplot2::element_text(size = legend_text),
      legend.key.size = ggplot2::unit(4,"mm"),
      legend.position = "top"
    )+
    ggplot2::theme(...) -> tn

  if(show_axis_title){
    tn =
      tn + ggplot2::theme(
      axis.title = element_text()
    )
  }
  tn
}


#' Using DT to print data.rame
#'
#' @param input_df input_data_frame
#' @param input_df_name data_frame_name
#' @param nrow how many rows each page
#' @param rownames whether showing the rownames
#' @param digits the digits of numeric column
#' @param fontSize font size, like "25%"
#'
#' @return DT obj
#' @export
#'
fun_utils_dt =
  function(input_df,input_df_name = "data_download",nrow = 5,rownames=F,digits = 3,fontSize = "25%"){
    num_column =    purrr::map_lgl(input_df, is.numeric)
    DT::formatRound(
      DT::formatStyle(
        DT::datatable(
          input_df,
          extensions = 'Buttons',
          width = "90%",
          rownames= rownames,
          options = list(
            paging = TRUE,
            searching = TRUE,
            pageLength = nrow,
            dom = 'Bfrtip',
            # extend = "csv",
            title ="data",
            scrollX = TRUE,
            buttons =list("copy",list(
              extend = 'csv',
              title = input_df_name))
          )),
        columns = colnames(input_df),
        fontSize = fontSize),
      columns = num_column,
      digits = digits
    )
  }



#' Conver tpval into star symbol
#'
#' @param input_pval the input value
#' @param ns non-significant labels
#'
#' @return a vector of strings
#' @export
#'
#' @examples
fun_utils_p2star = function(input_pval,ns = ""){
  purrr::map_chr(input_pval,function(x){
    dplyr::case_when(
      is.na(x)~NA_character_,
      x > 0.05~ns,
      x > 0.01~"*",
      x > 0.001~"**",
      x > 0.0001~"***",
      TRUE~"****"
    )
  })
}



#' eval the results
#'
#' @param input_texts the text inputs
#'
#' @return the val res
#' @export
#'
#' @examples
#' a = c("1+2", "1/2", "3")
#' fun_utils_eval(a)
fun_utils_eval = function(input_texts){
  purrr::map(input_texts,function(x){
    eval(parse(text = x))
  })

}


#' expression filter and/or quantitle normalization by row, if you want to do it by column, trannpose it and then transpose back
#'
#' @param input_mat input_mat
#' @param zero_pct zero_pct[0,1]
#' @param na_pct na_pct [0,1]
#' @param qn quantitle.normalzation by preprocess core
#'
#' @return filtered matrix
#' @export
#'
#' @examples
fun_utils_mat_filter = function(input_mat,zero_pct=0.2,na_pct = NULL,qn = F){
  if(!is.null(zero_pct)){
    # do zero pct
    zero_vaue = min(input_mat)
    minimun_count = ncol(input_mat) * zero_pct

    keep_index = rowSums(input_mat<=zero_vaue) < minimun_count

    input_mat = input_mat[keep_index,]
  }

  if(!is.null(na_pct)){

  }

  if(qn){
    input_mat = preprocessCore::normalize.quantiles(input_mat,keep.names = T)
  }

  input_mat
}




#' If the strins is to loner, say, 50 character, we would like to insert a new line into the strings, so that it could be printed beatifully
#'
#' @param input_strings  the stings list
#' @param string_check integer indicates the long strings
#' @param sep how the stings separate
#'
#' @return a list of stings with newline inserted
#' @export
#'
#' @examples
fun_utils_insert_newline_for_long_stings = function(
    input_strings,string_check = 50,sep = "_"){
  purrr::map_chr(input_strings,function(each_string){
    # 1) get sting length
    string_len = nchar(each_string)

    # 2) check the sting length
    if(string_len > string_check){
      string_start  = int(string_check/2)

      for(i in c(string_start : string_len)){
        if(stringr::str_sub(each_string,i,i) == sep){
          break
        }
      }

      paste0(stringr::str_sub(each_string,1,i-1),
             "\n",
             stringr::str_sub(each_string,i+1,string_len))

    }else{
      each_string
    }
  })

}



#' reinstall
#'
#' @param ...
#'
#' @export

fun_util_reinstall = function(...){
  detach("package:biodf")
  remove.packages("biodf")
  install.packages("~/data/project/pj/biodf_0.0.0.1.tar.gz",repos = NULL)
  .rs.restartR()
  require(biodf)
}


#' Get colors
#'
#' @param input_df input_df
#' @param input_variables input_variables
#' @param input_colors input_colors
#' @param color_shift color_shift
#'
#' @return a list with color map
#' @export
#'
#' @examples
fun_utils_get_color = function(
    input_df,input_variables=NULL,
    input_colors = NULL,
    color_shift = 4){

  # 1) init
  if(length(input_variables) == 0 ){
    input_variables = colnames(input_df)
  }
  if(length(input_colors) ==0){
    input_colors = c(
      ggpubr::get_palette("lancet",7),
      ggpubr::get_palette("jco",7),
      ggpubr::get_palette("default",7)
    )
  }else if(length(input_colors) == 1){
    input_colors = ggpubr::get_palette(input_colors,max(length(input_variables),10))
  }

  # 2) get color
  color_pointer = 0
  color_list = vector("list",length = length(input_variables))
  names(color_list) = input_variables

  for(i in seq_along(input_variables)){
    each_var = input_variables[i]
    # 1) get data
    var_data = na.omit(input_df[[each_var]])

    # 2) get color
    if(is.numeric(var_data)){
      # 1) get current color pointer
      color_pointer = color_pointer + 1
      color_tmp = input_colors[color_pointer]

      # 2) get numeric value
      min_val = min(var_data,na.rm = T)
      max_val = max(var_data,na.rm = T)

      color_low = grDevices::colorRampPalette(c("white",color_tmp))(10)[color_shift]
      color_list[[each_var]] =
        circlize::colorRamp2(
          breaks = c(min_val,max_val),
          colors = c(color_low,color_tmp))

    }else if(is.factor(var_data)){
      # 1) get color
      color_pointer = color_pointer + 1
      color_tmp = input_colors[color_pointer]

      # 2) get color gradient
      color_low = grDevices::colorRampPalette(c("white",color_tmp))(10)[color_shift]

      # 3) get data
      var_data = levels(var_data)

      # 4)
      color_list[[each_var]] =
        structure(grDevices::colorRampPalette(c(color_low,color_tmp))(length(var_data)),
                  names = var_data)


    }else{
      var_data = na.omit(unique(var_data))
      color_low_pointer = color_pointer+1
      color_high_pointer = color_pointer+ length(var_data)
      color_tmp = input_colors[c(color_low_pointer:color_high_pointer)]
      color_list[[each_var]] =
        structure(color_tmp,names = var_data)
      color_pointer = color_high_pointer
    }
  }
  color_list
}




#' Oncoprint data preparation
#'
#' @param input_df input_df
#' @param input_variables genes to visulize
#' @param input_variables_column gene column, like Hugo_Symbol in maf format
#' @param input_mut the mutation type column
#' @param input_sample the sample column, like Tumor_Sample_Barcode in maf format
#' @param input_sample_append input_sample_append
#' @param fun_aggre fun_aggre
#' @param input_col color for default mutation type
#' @param bg_col bg_col
#' @param cell_width cell_width
#' @param height_decrease height_decrease
#'
#' @return
#' @export
#'
#' @examples
fun_utils_df2oncoprint = function(
    input_df,
    input_variables,
    input_variables_column = "Hugo_Symbol",
    input_mut = "mutation_type",
    input_sample="sample",
    input_sample_append = NULL,
    fun_aggre = function(x) unique(x) %>% na.omit()%>% paste0(collapse = ";"),
    # alter function params
    input_col = NULL,
    bg_col ="gray",
    cell_height = 1,
    cell_width = 0.9,
    height_decrease = 0.05

){
  # 1) get dt
  input_dt =
    input_df %>%
    data.table::data.table() %>%
    data.table::dcast.data.table(
      formula = as.formula(paste0(input_sample,'~',input_variables_column)),
      value.var = input_mut,
      fun.aggregate = fun_aggre
    )

  # 2) append genes without mutation
  genes_not_found =input_variables %>%
    dplyr::setdiff(y = colnames(input_dt))

  for(i in seq_along(genes_not_found)){
    gene_tmp = genes_not_found[i]
    input_dt[[gene_tmp]] = ""
  }

  # 3) select the genes
  onco_mat =
    input_dt %>%
    dplyr::select(all_of(c(input_sample,input_variables))) %>%
    as.data.frame() %>%
    magrittr::set_rownames(.[[input_sample]]) %>%
    dplyr::select(-!!as.name(input_sample)) %>%
    t %>%
    as.data.frame()

  # 5) append samples
  input_sample_append =
    input_sample_append %>%
    dplyr::setdiff(colnames(onco_mat))

  for(i in seq_along(input_sample_append)){
    onco_mat[[input_sample_append[i]]] = ""
  }
  onco_mat = as.matrix(onco_mat)
  # 4) onco mut
  mut_type = input_df[[input_mut]]%>% na.omit() %>% unique
  if(is.null(input_col)){
    input_col = ggpubr::get_palette("Paired",length(mut_type))
  }else{
    stopifnot("input_col should be longer"= (length(input_col) >= length(mut_type)))
    input_col  = input_col[1:length(mut_type)]
  }
  # 5) get alter_type
  input_mut_type = c("background",mut_type)
  input_col = c(bg_col,input_col)
  purrr::map(seq_along(c(input_mut_type)),function(i){
    color_tmp = input_col[i]
    eval(rlang::expr(function(x, y, w, h) grid.rect(x, y, w*!!cell_width, h*(!!cell_height -!!height_decrease*!!i) ,gp = gpar(fill = !!color_tmp, col = NA))))
  }) %>%
    unlist(recursive = F) %>%
    purrr::set_names(input_mut_type) -> alter_list


  list(mat = onco_mat,
       alter_fun=alter_list,
       color_map = structure(input_col,names = input_mut_type))


}
