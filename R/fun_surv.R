#' Group numeric variables into high-low groups, based on the cutoff(Lower side included)
#'
#' @param input_df Input is data.frame-like object
#' @param input_variables Numeric variables
#' @param output_variables_names output names
#' @param input_pct range(0,1), If length(input_pct) == 1, using it as quantitle as cutoff;if length(input_pct) ==2, group higher than max percentage -> High, lower than the min percentage -> Low; If length(input_pct) >3, group data with minumu pvalue(survival, should have `time`,`status` columns)
#' @param label output labels, default: c("Low","High")
#' @param factor_level whether refactor group data, default is True
#'
#' @return A data.frame whose `input_variables` are grouped
#' @export
#'
#' @examples
#' data(tcga_clin)
#' clin_group = fun_surv_cutoff(
#'   tcga_clin,
#'   input_variables="age",
#'   output_variables_names = "Group",
#'   input_pct = 0.5)
fun_surv_cutoff = function(
    input_df,
    input_variables = "Risk_score",
    output_variables_names = "Group",
    input_pct=0.5,
    label = c("Low","High"),
    factor_level = T){
  if(length(input_pct) == 0){
    return(input_df)
  }else{
  stopifnot("`input_pct` should smaller than 1 and bigger than 0" = (min(input_pct) > 0 && max(input_pct) < 1))

  output_variables_names = fun_utils_broadcast(
    input_variables,
    output_variables_names)
  # 1) length == 1
  if(length(input_pct) == 1){
  for(i in seq_along(input_variables)){
    old_var = input_variables[i]
    new_var = output_variables_names[i]
    var_value = input_df[[old_var]]
    input_df[[new_var]] =
      ifelse(var_value <= quantile(var_value,input_pct,na.rm=T),
             label[1],
             label[2])
  }
  }else if(length(input_pct) == 2){
  # 2) length == 2
  for(i in seq_along(input_variables)){
    old_var = input_variables[i]
    new_var = output_variables_names[i]
    var_value = input_df[[old_var]]
    pct_min = min(input_pct)
    pct_max = max(input_pct)
    input_df[[new_var]] =
      dplyr::case_when(
        var_value <= quantile(var_value,pct_min,na.rm = FALSE)~label[1],
        var_value > quantile(var_value,pct_max,na.rm = FALSE)~label[2],
        TRUE~NA_character_)
  }
  }else {
  # 3) length > 3
  for(i in seq_along(input_variables)){
    old_var = input_variables[i]
    new_var = output_variables_names[i]
    var_value = input_df[[old_var]]

    pval_best_sorfar = 2
    pct_best_sofar = input_pct[1]

    for(pct in input_pct){
      df_tmp =fun_surv_cutoff(input_df,old_var,new_var,input_pct = pct)
      # cal the function again to geneate signle value
      pval_tmp = fun_surv_unicox(df_tmp,input_variables = new_var)
      pval_tmp = pval_tmp$pval
      pct_best_sofar =  ifelse(pval_tmp < pval_best_sorfar,pct,pct_best_sofar)
      pval_best_sorfar = ifelse(pval_tmp < pval_best_sorfar,pval_tmp,pval_best_sorfar)
    }

    message(paste0("Optimal cutoff point for ",old_var, " is: ",pct_best_sofar))
    input_df = fun_surv_cutoff(input_df,
                               old_var,
                               new_var,
                               input_pct  = pct_best_sofar,
                               label = label,
                               factor_level = factor_level)
  }
  }
  # 4) refactor
  if(factor_level){
    for(i in seq_along(output_variables_names)){
      var_new = output_variables_names[i]
      input_df[[var_new]] = factor(input_df[[var_new]],levels = label)
    }
  }

  input_df
  }

}

#' cutoff with quantile points, e.g. if you supplied input_pct as c(0.25,0.5,0.75), you will get four groups
#'
#' @param input_df  input_df
#' @param input_variables input_variables
#' @param output_variables_name output_variables_name
#' @param input_pct perentage should be in (0,1)
#' @param input_labels the lables length should be length(input_pct)+1
#' @param factor_level  T or F, whether to refactor as input_labels
#'
#' @return todo
#' @export
#'

fun_surv_cutoff_quantile = function(
    input_df,
    input_variables = "Risk_score",
    output_variables_name = "Group4",
    input_pct = c(0.25,0.5,0.75),
    input_labels = c("Low","Intermediate_low","Intermediate_high","High"),
    factor_level =T){
  # 1) get the input variables and output variables have the same length
  output_variables_name = fun_utils_broadcast(input_variables = input_variables,input_targets = output_variables_name)

  stopifnot("input_pct should be smaller than 1 and larger than 0" = (max(input_pct) <1 && min(input_pct) > 0))

  # 2) reorder the input_pct
  input_pct = sort(unique(input_pct))

  for(i in seq_along(input_variables)){
    # 1) get current variables
    current_var = input_variables[i]
    current_var_values = input_df[[current_var]]

    # 2) get current output variable names
    current_var_out_name = output_variables_name[i]
    input_df[[current_var_out_name]] = NA

    # 3) cutoff for each
    input_pct_values = quantile(current_var_values,input_pct,na.rm=T)
    input_pct_values = c(min(current_var_values,na.rm = T)-1,input_pct_values,max(current_var_values,na.rm=T)+1)

    # 4) for loop for the cutoff
    for(each_index in 1:(length(input_pct_values)-1)){
      current_pct_label = input_labels[each_index]
      current_low = input_pct_values[each_index]
      current_high = input_pct_values[(each_index+1)]

      input_df[[current_var_out_name]] =
        dplyr::case_when(
          (current_low < current_var_values & current_var_values <= current_high) ~ current_pct_label,
          TRUE~input_df[[current_var_out_name]]
        )
    }

    # 4) refactor
    if(factor_level){
      input_df[[current_var_out_name]] = factor(input_df[[current_var_out_name]],levels = input_labels)
    }
  }

  input_df
}


#' group var into three groups
#'
#'
#' @return df
#' @export

fun_surv_cutoff3 <- function(input_df,
                             input_variables = "Risk_score",
                             output_variables_name = "Group3",
                             input_pct = c(0.25,0.75),
                             input_labels = c("Low","Medium","High"),
                             factor_level =T){

  fun_surv_cutoff_quantile(input_df,
                           input_variables,
                           output_variables_name,
                           input_pct,
                           input_labels,
                           factor_level)

}


#' survival analysis
#'
#' @param input_df input data frame
#' @param input_variables input variables
#' @param palette color palette for KM
#' @param title KM plot showing title
#' @param risk_tables whether showing the risk tables
#' @param tables.theme the themes for the table
#' @param ggtheme the theme for all the
#' @param conf.int whether showing the the confident ribbon
#' @param xlab xlab
#' @param pval whether showing the pval
#' @param ... others params passing into the `ggsurvplot`
#'
#' @return ggsurvplot obj
#' @export
#'
#' @examples
#' data(tcga_clin)
#' fun_surv_km(head(tcga_clin,100),"Age_group")
fun_surv_km <- function(
    input_df,
    input_variables = "Group",
    input_pct = NULL,
    palette = "jco",
    title = NULL,
    risk_tables = T,
    tables.theme =ggplot2::theme(
      plot.title=ggplot2::element_blank()
    ),
    ggtheme = tn(),
    conf.int = F,
    xlab = "Months",
    pval=T,
    ...){

  # 1) map plot for each
  purrr::map(input_variables,function(x){
    # 1) rebuild the data
    data_tmp = data.frame(
      time = input_df$time,
      status = input_df$status,
      input_var = input_df[[x]]
    )
    data_tmp2 = fun_surv_cutoff(
      data_tmp,
      input_variables = "input_var",
      output_variables_names = "input_var",
      input_pct =input_pct
      )

    # 2) get the fit
    fit = survminer::surv_fit(
      formula  = survival::Surv(time,status)~input_var,
      data = data_tmp2
    )

    # 3) plot
    lenged.labs = sort(unique(data_tmp2$input_var))
    survminer::ggsurvplot(
      fit,
      surv.median.line = "hv",
      conf.int = conf.int,
      palette = palette,
      risk.table = risk_tables,
      tables.theme= tables.theme,
      legend.title = x,
      legend.labs =lenged.labs,
      title = title,
      xlab = xlab,
      pval=pval,
      ggtheme=ggtheme,
      ...
      ) -> p
    if(risk_tables){
    ggpubr::ggarrange(
      p$plot,p$table,
      nrow = 2,
      ncol=1,
      align = "v",
      heights = c(3,1.3))
    }else{
        p$plot
      }
  }) ->p_km_list

  if(length(p_km_list) == 1){
    p_km_list = p_km_list[[1]]
  }
  p_km_list
}

#' Parse survival time
#'
#' @param input_df input_df
#' @param input_surv_index three-length vector
#' @param input_na_value label some labels with NA value
#'
#' @return data.frame like obj with `time` and `status` column
#' @export
#'
fun_surv_parsing = function(
    input_df,
    input_surv_index=c("OS.time","OS","1"),
    input_na_value = NULL){

  if(!is.null(input_surv_index)){
    stopifnot("`input_surv_index` should be `NULL` or a vector with three elements\n for example c(OS.time,status,dead)"=(length(input_surv_index) == 3))
    time_tmp = input_surv_index[1]
    status_tmp = input_surv_index[2]
    status_l = input_surv_index[3]
    # replace na label with NA value in R
    if(!is.null(input_na_value)){
      input_df[[status_tmp]] = ifelse(input_df[[status_tmp]] %in% input_na_value,
                                      NA_character_,
                                      input_df[[status_tmp]])

    }
    stopifnot("Status with two unique labels are considered normal" =
                length(na.omit(unique(input_df[[status_tmp]]))) == 2)
    input_df$time = input_df[[time_tmp]]
    input_df$status = input_df[[status_tmp]]
    input_df$status = ifelse(input_df$status == status_l,1,0)
  }
  return(input_df)
}

#' Unicox for survival model, make sure you have time and status column, and label properly
#'
#' @param input_df Input data.frame, with variables and time, status columns
#' @param input_variables Input variables
#'
#' @return data.frame-like objects
#' @export
#'
#' @examples
#' data(tcga_clin)
#' fun_surv_unicox(head(tcga_clin,100),"Age_group")
#'
fun_surv_unicox = function(
    input_df,
    input_variables = "Risk_score",
    input_pct = NULL
){

  purrr::map_df(input_variables,function(x){
    if(is.numeric(input_df[[x]])){
      input_df = fun_surv_cutoff(input_df,x,x,
                                 input_pct = input_pct,
                                 label = c("Low","High"),
                                 factor_level = T)
    }
    # 1) formula
    data_tmp = data.frame(
      time = input_df$time,
      status = input_df$status,
      input_var = input_df[[x]]
    )

    # 2) analysis
    res = survival::coxph(
      formula=survival::Surv(time,status)~input_var,
      data = data_tmp)


    # 3) add logrank pvalue
    if(!is.numeric(data_tmp$input_var)){
      surv_fit =
        survminer::surv_fit(
          survival::Surv(time,status)~input_var,
          data =data_tmp
        )
      km_res = survminer::surv_pvalue(surv_fit)

    }else{
      km_res = list()
      km_res$pval = NA
    }


    # 3) return df
    df1 = as.data.frame(coef(summary(res)))
    rownames(df1) = stringr::str_replace(rownames(df1),
                                         "^input_var",
                                         x)
    df1$features = rownames(df1)
    df1$hr = df1$`exp(coef)`
    df1$ci_low = exp(df1$coef - 1.96 * df1$`se(coef)`)
    df1$ci_up= exp(df1$coef + 1.96 * df1$`se(coef)`)
    df1$pval = df1[[5]]
    df1$pval_logrank = NA
    df1$pval_logrank[1] = km_res$pval
    df1 = df1[,c("features","hr","ci_low","ci_up","pval","pval_logrank"),drop =F]
    if(length(input_pct) != 0){
      rownames(df1) = stringr::str_remove(rownames(df1),"High$")
      df1$features = stringr::str_remove(df1$features,"High$")
    }
    df1
  })
}

#' Multicox survival model, make sure you have time and status column, and label properly
#'
#' @param input_df input data.frame, with variables and time, status columns
#' @param input_variables input variables
#'
#' @return data.frame-like objects
#' @export
#'
#' @examples
#' data(tcga_clin)
#' fun_surv_multicox(head(tcga_clin,100),c("age","Age_group"))
#
fun_surv_multicox = function(
    input_df,
    input_variables = "Risk_score"){

  # 1)get all the data
  input_df_tmp <-
    dplyr::select(input_df,all_of(c(input_variables,"time","status")))

  # 2) get the fit
  fit = survival::Surv(time,status)~.

  # 3) cox plot
  res = survival::coxph(fit,data = input_df_tmp)

  # 4) change the results into data.frame format
  df1 = as.data.frame(coef(summary(res)))
  rownames(df1) = stringr::str_remove_all(rownames(df1),"`")

  # 5)
  df1$features = rownames(df1)
  df1$hr = df1$`exp(coef)`
  df1$ci_low = exp(df1$coef - 1.96 * df1$`se(coef)`)
  df1$ci_up= exp(df1$coef + 1.96 * df1$`se(coef)`)
  df1$pval = df1[[5]]
  df1[,c("features","hr","ci_low","ci_up","pval"),drop =F]
}

#' roc analysis for survival data
#'
#' @param input_df input.frame data
#' @param input_variables input variables
#' @param input_times times for survival time analysis
#' @param input_ds_name ds_names
#' @param merge_roc_auc merge input variables roc/auc in one data.frame respectively
#'
#' @return if `merge_roc_auc` is `TRUE`,1)list with two data.frame(roc and auc) 2)Otherwise, lists with same length to input_variables
#' @export
#'
#' @examples
#' # merge the output into two dataframe
#' res_list =
#'  fun_surv_roc(
#'  tcga_clin,
#'  input_variables = "age",
#'  merge_roc_auc=T )
#'

fun_surv_roc =function(
    input_df,
    input_variables = "Risk_score",
    input_times = c(1,2,3,4,5)*12,
    input_ds_name = "Not supplied",
    merge_roc_auc =T){
  suppressMessages({
    require(survival)
    require(timeROC)
  })
  roc_list =
    purrr::map(input_variables,function(x){
    # 1) fit timeROC
    surv_roc = timeROC::timeROC(T = input_df$time,
                           delta = input_df$status,
                           marker = input_df[[x]],
                           cause = 1,
                           times = input_times,
                           iid = T)
    # 2) get roc data.frmae
    # timeROC results parsing
    dat = data.frame(fpr = as.numeric(surv_roc$FP),
                     tpr = as.numeric(surv_roc$TP),
                     time = rep(as.factor(input_times),
                                each = nrow(surv_roc$TP)))
    dat$features =x
    dat$ds_name = input_ds_name

    # 3) get auc data.frame
    auc_df = data.frame(time = input_times,
                        auc = surv_roc$AUC)

    # add confint
    auc_ci = confint(surv_roc)[[1]]
    auc_df$ci_low = c(as.numeric(auc_ci[,1])/100)[1:length(input_times)]
    auc_df$ci_up = c(as.numeric(auc_ci[,2])/100)[1:length(input_times)]
    auc_df$features = x
    auc_df$ds_name = input_ds_name
    rownames(auc_df) = paste0(x,"_",rownames(auc_df))
    return(list(roc = dat,
                auc = auc_df))
  })

  # 2) set names
  names(roc_list) = input_variables

  # 3) return
  if(merge_roc_auc){
    # 1) retur merged data.frame
    roc_df = purrr::map_df(roc_list,function(x){
      x$roc
    })
    auc_df = purrr::map_df(roc_list,function(x){
      x$auc
    })
    list(roc = roc_df,
         auc = auc_df)
  }else{
    #2) return raw data
    roc_list
  }
}


#' Plot time roc by the results of `fun_surv_roc`, extract the roc and auc
#'
#' @param input_surv_roc_list input_surv_roc_list
#' @param palette using get_paltte function in ggpubr, if its length equals to 1
#' @param title title
#' @param legend_suffix legend_suffix
#' @param legend_prefix legend_prefix
#' @param legend.position legend.position
#'
#' @return ggplot list
#' @export
#'

fun_plot_surv_roc_time = function(
    input_surv_roc_list,
    palette = "jco",
    title = "",
    legend_suffix = "-m survival ",
    legend_prefix = "AUC of ",
    legend.position = c(0.6,0.2)){
  dat  = input_surv_roc_list$roc
  auc = input_surv_roc_list$auc

  # 1) palette parsing
  if(length(palette) == 1){
    time_n = nrow(auc)
    palette = ggpubr::get_palette(palette = palette,k = time_n)
  }

  p =
    ggplot2::ggplot(data = dat) +
    ggplot2::geom_abline(slope = 1,color = "gray") +
    ggplot2::geom_line(aes(x= fpr,y=tpr,color=time))+
    # legend
    ggplot2::scale_color_manual(
      name = NULL,values = palette,
      labels = paste0(legend_prefix,
                      auc$time,
                      legend_suffix,
                      format(round(auc$auc,2),nsmall = 2)))+
    # line
    ggplot2::theme_bw()+
    ggplot2::theme(panel.grid =ggplot2::element_blank(),
                   legend.background = ggplot2::element_blank(),
                   legend.position = legend.position )+
    ggplot2::labs(
      x = "1 - Specificity",
      y = "Sensitivity",
      title=title)+
    ggplot2::scale_x_continuous(expand = c(0.005,0.005))+
    ggplot2::scale_y_continuous(expand = c(0.005,0.005))+
    ggplot2::coord_fixed()

  return(p)
}

ft_surv_concordance_fit_to_df = function(input_fit,fit_now){
  stats_tmp = abs(qt(p = 0.025,df = input_fit$n -1 ))
  conf_value = stats_tmp * sqrt(input_fit$var)
  data.frame(row.names = paste0(fit_now,collapse = "&"),
             concordance = input_fit$concordance,
             ci_low = input_fit$concordance -  conf_value,
             ci_up = input_fit$concordance + conf_value)
}

ft_concordance_fit_pval = function(input_fit1,input_fit2){
  m_diff = input_fit1$concordance - input_fit2$concordance
  sd_new = sqrt(input_fit1$var/input_fit1$n +input_fit2$var/input_fit2$n)
  t = -abs(m_diff/sd_new)
  pt(q = t,df = 1)

}

#' Calculate concordance step
#'
#' @param input_df input_df
#' @param input_variables input_variables
#' @param input_ref set 1
#'
#' @return a data.frame with value
#' @export
#'

fun_surv_concordance = function(input_df,
                                input_variables,
                                input_ref = as.integer(1)){
  warnings("Pvalue calculation still need to be done")
  stopifnot("input_ref must be integer index"=is.integer(input_ref))

  # 1) get reference fit
  fit_list =
    purrr::map(input_variables,function(x){
      df_tmp = dplyr::select(input_df,time,status,dplyr::all_of(x))
      survival::concordance(survival::coxph(survival::Surv(time,status)~.,data = df_tmp))
    })

  # 2)
  fit_ref = fit_list[[input_ref]]
  df_init = ft_surv_concordance_fit_to_df(fit_ref,input_variables[[input_ref]])
  df_init$pval = NA
  for(i in seq_along(fit_list)){
    if(i != input_ref){
      df_new =  ft_surv_concordance_fit_to_df(fit_list[[i]],input_variables[[i]])
      df_new$pval = ft_concordance_fit_pval(fit_ref,fit_list[[i]])
      df_init = rbind(df_init,df_new)

    }
  }
  df_init
}

#' same as fun_surv_km
#'
#' @param ... passing to fun_surv_km
#'
#' @export
#'
fun_plot_surv_km <- function(...){
  input_args = list(...)
  do.call(fun_surv_km,input_args)
}



#' Plot unicox pvalues and HR with differnet cutoff
#'
#' @param input_df input_df
#' @param input_variables input_variables
#' @param input_pct input_pct,deafulat from 0.2 to 0.8 by 0.01
#' @param axis_fold axis_fold is 2
#' @param x_label x_label
#' @param y_label y_label
#' @param sec_y_axis_label sec_y_axis_label
#' @param grid_color grid_color
#' @param vline_lty vline_lty
#' @param vline_color vline_color
#' @param line_color line_color
#'
#' @return todo
#' @export
#'
fun_plot_surv_cutoff <- function(input_df,
                                 input_variables = "Risk_score",
                                 input_pct = seq(0.2,0.8,0.01),
                                 axis_fold = 2,
                                 return_cutoff = F,
                                 x_label = "Quantile(%)",
                                 y_label = "P-value",
                                 sec_y_axis_label = "Hazard ratio",
                                 grid_color = "gray90",
                                 vline_lty = "twodash",
                                 vline_color = "#FF9B82",
                                 line_color = c("#8BE8E5","#9400FF")
                                 ){

  df <-
    purrr::map_df(seq(0.1,0.9,0.1),\(each_pct){
      input_df %>%
        fun_surv_cutoff(input_variables,input_variables,input_pct = each_pct) %>%
        fun_surv_unicox(input_variables) %>%
        dplyr::mutate(Quantile = each_pct * 100)
    }) %>%
    dplyr::select(Quantile,hr,pval) %>%
    dplyr::mutate(hr = hr / axis_fold) %>%
    tidyr::pivot_longer(-Quantile)

  xintercept = df %>%
    dplyr::filter(name == "pval") %>%
    dplyr::filter(value == min(value)) %>%
    dplyr::pull(Quantile)

  df%>%
    ggplot2::ggplot(ggplot2::aes(Quantile,value,group = name))+
    ggplot2::geom_line(ggplot2::aes(color = name,lty = name))+
    ggplot2::geom_point(ggplot2::aes(shape = name))+
    tn()+
    ggplot2::scale_y_continuous(
      sec.axis = ggplot2::sec_axis(~.* axis_fold,name = sec_y_axis_label))+
    ggplot2::theme(panel.grid.major = ggplot2::element_line(color =grid_color,linetype = "dashed"))+
    ggplot2::geom_vline(xintercept = xintercept,lty = vline_lty,color = vline_color)+
    ggplot2::scale_color_manual(values = line_color)+
    ggplot2::labs(x = x_label, y = y_label) -> p
  if(return_cutoff){
    list(p = p, cutoff = xintercept)
  }else{
    p
  }
}



#' clinical data survival parsing quick entry
#'
#' @param input_df clinical data
#' @param shortcut which data to parse
#' @param fun_surv_parsing_arge args passing to fun_surv_parsing
#'
#' @return data.frame: clinical data
#' @export
#'
fun_surv_parsing_shortcut <- function(
    input_df,
    shortcut = 'pfs',
    fun_surv_parsing_arge =list()){
  fun_surv_parsing_arge$input_df = input_df
  # 1) get shortcuts
  shortcut =match.arg(shortcut,c('pfs','dss','pfi','dfi',
                                 'dfs',
                                 'os_cbio','dfs_cbio',
                                 'os_Cbio','dfs_Cbio'))

  # 2) get the surv index
  fun_surv_parsing_arge$input_surv_index =
    switch (shortcut,
      pfs = c('pfs_time','pfs_status',1),
      dss = c('dss_time','dss_status',1),
      pfi = c('pfi_time','pfi_status',1),
      dfi = c('dfi_time','dfi_status',1),
      dfs = c('dfs_time','dfs_status',1),

      os_cbio = c('overall_survival_months','overall_survival_status','1:DECEASED'),
      dfs_cbio = c('disease_free_months','disease_free_status','1:Recurred/Progressed'),

      os_Cbio = c('Overall Survival (Months)','Overall Survival Status','1:DECEASED'),
      dfs_Cbio = c('Disease Free (Months)','Disease Free Status','1:Recurred/Progressed'),

    )
  # results
  do.call(fun_surv_parsing,fun_surv_parsing_arge)
}



#' Nomogram
#'
#' @param input_df input_df
#' @param input_variables input_variables
#' @param times times
#' @param B boostrap
#' @param m m
#' @param lim lim
#' @param time_suffix months
#' @param outfile pdf save postion
#'
#' @return todo
#' @export
#'
fun_surv_nomogram <-  function(input_df,input_variables,times = c(12,24,36,48,60), B=10,m=50,lim = c(0,1),time_suffix = "-month",outfile = "./plot/nomogram"){
  suppressMessages(library(riskRegression))
  suppressMessages(library(regplot))
  suppressMessages(library(dplyr))
  suppressMessages(library(magrittr))
  suppressMessages(library(rms))
  suppressMessages(library(survival))
  suppressMessages(library(survcomp))

  # 1) get data
  input_df %<>%
    dplyr::select(dplyr::all_of(c("time","status",input_variables))) %>%
    data.frame(check.names = T) %>%
    #filter(time > 0) %>%
    na.omit() %>%
    {.}
  input_variables = colnames(input_df) %>%
    dplyr::setdiff(c("time","status"))
  # 2) lasso_cox
  formular_coxph = as.formula(paste0(
    "Surv(time,status)~",paste0(colnames(input_df)[-c(1,2)],collapse = "+")
  ))
  cox_model = coxph(formular_coxph,data = input_df,x=T,y=T)
  # 3) nomogram
  regplot(cox_model,
          failtime = times,
          prfail = F,
          title = "Nomogram",
          points = T,
          center = F,
          subticks = F) -> nomo_results

  # 4) concordance
  #concordance(cox_model,data = input_df) %>% print
  #nomo.output=predict(cox_model)
  cindex <- concordance.index(predict(cox_model),
                              surv.time = input_df$time,
                              surv.event = input_df$status,method = "noether")
  paste0("C_index\n",
         round(cindex$c.index,2),"\n",
         "C_index_confint\n",
         round(cindex$lower,2),"-",
         round(cindex$upper,2)) %>% message()
  #
  # # 4)prediction
  input_df[["Points"]] = 0
  for(f in seq_along(input_variables)){
    df = nomo_results[[f]] %>% as.data.frame()
    feature_now = intersect(input_variables,names(nomo_results[[f]]))
    # numeric lm
    if(is.numeric(df[[feature_now]])){
      for_mula = as.formula(
        paste0("Points~",feature_now)
      )
      reg = glm(for_mula,data = df)
      newdta = data.frame(time = input_df[[1]])
      newdta[[feature_now]] = input_df[[feature_now]]
      newdta[[feature_now]] = as.numeric(newdta[[feature_now]])
      rest =  predict(reg,newdata = newdta)
    }else{
      feature_data =df[[feature_now]]
      feature_point = df[["Points"]]
      feature_point[is.na(feature_point)] = 0
      newdta = data.frame(time = input_df[[1]])
      newdta[[feature_now]] = input_df[[feature_now]]
      # map data with points
      newdta$points = data.frame(row.names = feature_data,
                                 score = feature_point)[newdta[[feature_now]],1]
      rest = as.numeric(newdta$points)
    }
    #
    input_df[["Points"]] = input_df[["Points"]]  + rest
  }

  # 5) timeroc
  input_df %>%
    fun_surv_roc(input_variables = "Points",
                 input_times =  times) %>%
    fun_plot_surv_roc_time() -> p_roc
  par(mar = c(5,4,3,2))
  par(oma = c(0,0,0,0))
  #par(plt = c(5,4,4,2))

  # 6) calibration curve
  pdf(paste0(outfile,"_cindex.pdf"))

  for(time in times){
    set.seed(1)
    print(time)
    f <- cph(formular_coxph, x=T, y=T, surv=T,data=input_df, time.inc=time)
    #print(f)
    cal = calibrate(f, u=time,cmethod="KM",method = "boot", m=m, B=B)
    #print(cal)
    plot(cal,xlim=lim,ylim=lim,xlab="Predicted survival",ylab="Actual survival",subtitles=F)
    title(paste0(time,time_suffix," calibration curve"),cex.main = 0.7)
  }
  dev.off()


  return(list(
    data = input_df,
    time_roc = p_roc,
    fit = cox_model
  ))
}
