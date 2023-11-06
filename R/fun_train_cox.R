#' Cox training
#'
#' @param input_df input_ff
#' @param input_variables variables for cox regression
#' @param input_y time and status column
#' @param alpha 1->lasso, 0 -> ridge elastic (0,1)
#' @param seed reproducible seed
#' @param return_fit whether return fit
#' @param lambda coefficients lambda
#'
#' @return coefficients coef or list with fit
#' @export
#'
fun_train_cox = function(input_df,
                         input_variables,
                         input_y = c("time","status"),
                         alpha = 1,
                         return_fit = F,
                         seed = 1,
                         lambda="lambda.min"){
  input_df = fun_surv_data_select = fun_surv_select(input_df,input_variables ,input_y)

  x = as.matrix(dplyr::select(input_df,dplyr::all_of(input_variables)))
  y = as.matrix(dplyr::select(input_df,all_of(input_y)))
  set.seed(seed)
  fit = glmnet::cv.glmnet(x = x,y = y,alpha = alpha,family = "cox")
  coef_data = glmnet::coef.glmnet(fit,s=lambda)
  variable_index = coef_data@i + 1
  variables_final = coef_data@Dimnames[[1]][variable_index]
  variables_value = coef_data@x

  coef_final =
    data.frame(row.names = variables_final,
               coef = variables_value)

  if(return_fit){
    list(fit = fit,
         coef = coef_final)
  }else{
    coef_final
  }
}

#' Step cox analysis
#'
#' @param input_df input_df
#' @param input_variables input_variables
#' @param direction `both`, `forward` or `backward`
#' @param input_y `time` and `status` column, please set to `time` and `status` all the time
#' @param return_fit whether return fit
#'
#' @return a list with fit and coef, or just coef only
#' @export
#'
fun_train_cox_step = function(input_df,
                              input_variables,
                              direction = "forward",
                              input_y = c("time","status"),
                              return_fit = F){
  input_df = fun_surv_data_select = fun_surv_select(input_df,input_variables ,input_y)


  # 2) get formula
  surv_formula = survival::Surv(time,status)~.

  # 3) get  results
  cox_results_multi = survival::coxph(surv_formula,input_df)

  # 4) step
  step_cox = MASS::stepAIC(cox_results_multi,
                           direction = direction,
                           trace=F)

  # 5) coefficients
  coef_final = as.data.frame(coef(step_cox)) %>%
    magrittr::set_rownames(rownames(.),str_remove_all(rownames(.),"`"))
  colnames(coef_final) = "coef"

  if(return_fit){
    list(fit = step_cox,
         coef = coef_final)
  }else{
    coef_final
  }
}

#' Survival random forest
#'
#' @param input_df input_df
#' @param input_variables input_variables
#' @param input_y `time`,`status`
#' @param nsplit pass to `randomForestSRC::rfsrc`
#' @param nodesize pass to `randomForestSRC::rfsrc`
#' @param conservative pass to `randomForestSRC::rfsrc`
#' @param return_fit whether return fit
#' @param seed reproducible seed
#'
#' @return filtered variables or list with fit and filtered variables
#' @export
#'
fun_train_cox_rf = function(input_df,
                            input_variables,
                            input_y = c("time","status"),
                            nsplit = 50,
                            nodesize = 5,
                            conservative = "low",
                            return_fit = F,
                            seed = 1){

  # 1) get data
  input_df = fun_surv_data_select = fun_surv_select(input_df,input_variables ,input_y)


  # 3)
  Surv = survival::Surv
  sur_for = Surv(time,status)~.
  obj = randomForestSRC::rfsrc(
    sur_for,
    input_df,
    seed = seed,
    ntree = 1000,
    nodesize = nodesize,
    nsplit = nsplit,
    importance =TRUE)
  vs.pbc <- randomForestSRC::var.select(object = obj,
                                        conservative = conservative)
  topvars = as.character(vs.pbc$topvars)
  if(return_fit){
    list(fit = obj,
         topvars = topvars)
  }else{
    topvars
  }

}


#' CoxBoost
#'
#' @param input_df input_df
#' @param input_variables input_variables
#' @param input_y `time` and `status`
#' @param return_fit whether return fit
#' @param seed reproducible seed
#'
#' @return coef aor list with fit and coef
#' @export
#'
fun_train_cox_boost = function(input_df,
                               input_variables,
                               input_y = c("time","status"),
                               return_fit = F,
                               seed = 1){
  suppressMessages(require(CoxBoost))

  # 1) get data
  input_df = fun_surv_data_select = fun_surv_select(input_df,input_variables ,input_y)

  # get data
  time = input_df[[input_y[1]]]
  status = input_df[[input_y[2]]]
  x = as.matrix(dplyr::select(input_df,dplyr::one_of(input_variables)))
  # 1) cox penalty
  set.seed(seed)
  cox_penalty = CoxBoost::optimCoxBoostPenalty(
    time=time,
    status=status,
    x = x,
    trace=F,
    start.penalty = 500
  )
  set.seed(seed)
  cv_cox_step = CoxBoost::cv.CoxBoost(time = time,
                                      status = status,
                                      x = x,
                                      maxstepno=cox_penalty$cv.res$optimal.step,
                                      K=10,
                                      type="verweij",
                                      penalty=cox_penalty$penalty
  )
  set.seed(seed)
  cb_fit <-CoxBoost::CoxBoost(time = time,
                              status =status,
                              x =x,
                              stepno=cv_cox_step$optimal.step,
                              penalty=cox_penalty$penalty)

  coef_final = as.data.frame(coef(cb_fit)) %>%
    magrittr::set_rownames(rownames(.),str_remove_all(rownames(.),"`"))
  colnames(coef_final) = "coef"
  coef_final = dplyr::filter(coef_final,coef!=0)

  # return
  if(return_fit){
    list(fit = cb_fit,
         coef = coef_final)
  }else{
    coef_final

}}



fun_surv_select = function(input_df,input_variables,input_y=c("time","status")){
  input_df = dplyr::select(input_df,dplyr::all_of(c(input_variables,input_y)))
  input_df = input_df[input_df[[input_y[1]]] > 0,]
  input_df = na.omit(input_df)
  input_df
}




#' Survival feature selection
#'
#' @param input_df input_df
#' @param input_variables input_variables
#' @param outfile the results saved on disk
#' @param pval_cutoff pval_cutoff for unicox and multicox
#' @param methods_remove not using the methods
#' @param methods "lasso","eNet75","eNet50","eNet25","stepForward","stepBackward","stepBoth","multicox","unicox","boost","rf"
#' @param fun_train_cox_args fun_train_cox_args
#' @param fun_train_cox_step_args fun_train_cox_step_args
#' @param fun_surv_multicox_args fun_surv_multicox_args
#' @param fun_surv_unicox_args fun_surv_unicox_args
#' @param fun_train_cox_boost_args fun_train_cox_boost_args
#' @param fun_train_cox_rf_args fun_train_cox_rf_args
#'
#' @return feature list
#' @export
#'
#' @examples
fun_feature_selection_surv <- function(
    input_df,input_variables,
    outfile,
    pval_cutoff = 0.05,
    methods_remove = NULL,
    methods = c(
      "lasso","eNet75","eNet50","eNet25",
      "stepForward","stepBackward","stepBoth",
      "multicox","unicox","boost","rf"),
    # args
    fun_train_cox_args = list(),
    fun_train_cox_step_args = list(),
    fun_surv_multicox_args = list(),
    fun_surv_unicox_args = list(),
    fun_train_cox_boost_args = list(),
    fun_train_cox_rf_args = list()
){

  if(!file.exists(outfile)){
    message("results not found, running...")

    # 1) get data
    df_used <-
      input_df %>%
      dplyr::select(all_of(input_variables),time,status) %>%
      dplyr::filter(time>0)

    methods_select <- methods %>%
      dplyr::setdiff(methods_remove)

    message("Current select feature selection methods: ",
            length(methods_select),
            "\n",
            paste0(methods_select,collapse = ","),
            "\nCurrent select feature selection features:",length(input_variables))

    # args parsing
    fun_train_cox_args$input_df = df_used
    fun_train_cox_args$input_variables = input_variables

    fun_train_cox_step_args$input_df = df_used
    fun_train_cox_step_args$input_variables = input_variables

    fun_surv_multicox_args$input_df = df_used
    fun_surv_multicox_args$input_variables = input_variables

    fun_surv_unicox_args$input_df = df_used
    fun_surv_unicox_args$input_variables = input_variables

    fun_train_cox_boost_args$input_df = df_used
    fun_train_cox_boost_args$input_variables = input_variables

    fun_train_cox_rf_args$input_df = df_used
    fun_train_cox_rf_args$input_variables = input_variables

    # 2) methods train step
    tictoc::tic()
    results <- purrr::map(methods_select,function(each_method){
      message("-------><-------")
      message(each_method,": runing...")
      if(each_method == "lasso"){
        fun_train_cox_args$alpha = 1
        features_remained <- do.call(fun_train_cox,fun_train_cox_args) %>% rownames()
      }else if(each_method == "eNet75"){
        fun_train_cox_args$alpha = 0.75
        features_remained <- do.call(fun_train_cox,fun_train_cox_args) %>% rownames()
      }else if(each_method == "eNet50"){
        fun_train_cox_args$alpha = 0.5
        features_remained <- do.call(fun_train_cox,fun_train_cox_args) %>% rownames()
      }else if(each_method == "eNet25"){
        fun_train_cox_args$alpha = 0.25
        features_remained <- do.call(fun_train_cox,fun_train_cox_args) %>% rownames()
      }else if(each_method == "stepForward"){
        fun_train_cox_step_args$direction = "forward"
        features_remained <- do.call(fun_train_cox_step,fun_train_cox_step_args) %>% rownames()
      }else if(each_method == "stepBackward"){
        fun_train_cox_step_args$direction = "backward"
        features_remained <- do.call(fun_train_cox_step,fun_train_cox_step_args) %>% rownames()
      }else if(each_method == "stepBoth"){
        fun_train_cox_step_args$direction = "both"
        features_remained <- do.call(fun_train_cox_step,fun_train_cox_step_args) %>% rownames()
      }else if(each_method == "multicox"){
        message("cutoff: ", pval_cutoff)
        features_remained <- do.call(fun_surv_multicox,fun_surv_multicox_args) %>%
          dplyr::filter(pval < pval_cutoff) %>%
          rownames()
      }else if(each_method == "unicox"){
        message("cutoff: ", pval_cutoff)
        features_remained <- do.call(fun_surv_unicox,fun_surv_unicox_args) %>%
          dplyr::filter(pval < pval_cutoff) %>%
          rownames()
      }else if(each_method == "boost"){
        features_remained <- do.call(fun_train_cox_boost,fun_train_cox_boost_args) %>% rownames()
      }else if(each_method == "rf"){
        features_remained <- do.call(fun_train_cox_rf,fun_train_cox_rf_args)
      }else{
        warning("invalida input: continue next:")
        warning("This is an invalid method input,check the input,or remove this methods")
        features_remained = "This is an invalid method input,check the input,or remove this methods"
      }
      # print final kept feature
      message(each_method,"kept features: ",length(features_remained))
      return(features_remained)
    }) %>%
      set_names(methods_select)
    tictoc::toc()
    results %>% saveRDS(outfile)
  }
  message("Loading results")
  results = readRDS(outfile)
  for(each_method in names(results)){
    message(each_method,": ",length(results[[each_method]]))
  }
  return(results)
}


#' Train step, use with fun_feature_selection_surv and fun_design_train_combn
#'
#' @param input_df input_df
#' @param design_df two column df, feature selection methods and training methods
#' @param input_features_list the feature selected methods list
#' @param outfile the results saving outfile
#' @param feature_selection_method the colname of feature selection in  design_df
#' @param train_method the colname of feature selection in  train method
#' @param fun_train_cox_args fun_train_cox
#' @param fun_train_cox_step_args fun_train_cox_args
#' @param fun_surv_multicox_args fun_surv_multicox
#' @param fun_train_cox_boost_args fun_train_cox_boost
#'
#' @return list of model coef
#' @export
#'
#' @examples
fun_train_multiple_cox <- \(input_df,
                        design_df,
                        input_features_list,
                        outfile,
                        feature_selection_method = "m1",
                        train_method = "m2",
                        # args
                        fun_train_cox_args = list(),
                        fun_train_cox_step_args = list(),
                        fun_surv_multicox_args = list(),
                        fun_train_cox_boost_args = list()
){
  if(!file.exists(outfile)){
    message("model not found, running...")
    results = list()
    # remoing all ticks
    input_features_list =
      purrr::map(input_features_list,\(each_var){
        each_var %>% stringr::str_remove_all("`")
      })
    all_var_check = input_features_list %>%
      unlist() %>%
      unname() %>%
      unique()
    df_check <-  input_df %>%
      dplyr::select(dplyr::all_of(all_var_check))
    fun_train_cox_args$input_df = input_df
    fun_train_cox_step_args$input_df = input_df
    fun_surv_multicox_args$input_df = input_df
    fun_train_cox_boost_args$input_df = input_df

    for(i in seq_along(design_df[[1]])){
      m1 = design_df[[feature_selection_method]][i]
      each_method = design_df[[train_method]][i]
      model_name = paste0(m1,"&",each_method)
      message("-------------><-------------")
      message("model_name: ",model_name)
      # get genes
      input_variables = input_features_list[[m1]]
      # get args
      # args parsing
      fun_train_cox_args$input_variables = input_variables
      fun_train_cox_step_args$input_variables = input_variables
      fun_surv_multicox_args$input_variables = input_variables
      fun_train_cox_boost_args$input_variables = input_variables

      # train step
      if(each_method == "lasso"){
        fun_train_cox_args$alpha = 1
        model_df <- do.call(fun_train_cox,fun_train_cox_args)
      }else if(each_method == "eNet75"){
        fun_train_cox_args$alpha = 0.75
        model_df <- do.call(fun_train_cox,fun_train_cox_args)
      }else if(each_method == "eNet50"){
        fun_train_cox_args$alpha = 0.5
        model_df <- do.call(fun_train_cox,fun_train_cox_args)
      }else if(each_method == "eNet25"){
        fun_train_cox_args$alpha = 0.25
        model_df <- do.call(fun_train_cox,fun_train_cox_args)
      }else if(each_method == "ridge"){
        fun_train_cox_args$alpha = 0
        model_df <- do.call(fun_train_cox,fun_train_cox_args)
      }else if(each_method == "stepForward"){
        fun_train_cox_step_args$direction = "forward"
        model_df <- do.call(fun_train_cox_step,fun_train_cox_step_args)
      }else if(each_method == "stepBackward"){
        fun_train_cox_step_args$direction = "backward"
        model_df <- do.call(fun_train_cox_step,fun_train_cox_step_args)
      }else if(each_method == "stepBoth"){
        fun_train_cox_step_args$direction = "both"
        model_df <- do.call(fun_train_cox_step,fun_train_cox_step_args)
      }else if(each_method == "multicox"){
        model_df <- do.call(fun_surv_multicox,fun_surv_multicox_args) %>%
          as.data.frame() %>%
          magrittr::set_rownames(.$features) %>%
          dplyr::mutate(coef = log(hr)) %>%
          dplyr::select(coef)
      }else if(each_method == "boost"){
        model_df <- do.call(fun_train_cox_boost,fun_train_cox_boost_args)
      }else {
        warning("invalida input: continue next:")
        warning("This is an invalid method input,check the input,or remove this methods")
        model_df = "This is an invalid method input,check the input,or remove this methods"
      }
      message("model variables count",": ",nrow(model_df))
      message("model variables: ",paste0(rownames(model_df),collapse = ","))

      results[[model_name]] = model_df
    }
    results %>% saveRDS(outfile)
  }
  results = readRDS(outfile)
  for(model_name in names(results)){
    model_df = results[[model_name]]
    message("model_name",": ",model_name)
    message("model variables count",": ",nrow(model_df))
    message("model variables: ",paste0(rownames(model_df),collapse = ","))
  }
  results
}

#' Feature selection methods
#'
#' @param input_methods cancidate methods
#' @param notUsedForTraining not used as train
#'
#' @return df
#' @export
#'
#' @examples
fun_design_cox <- function(input_methods,notUsedForTraining = c("rf","unicox")){
  purrr::map(input_methods,\(each_name){
    purrr::map(input_methods %>% setdiff(notUsedForTraining),function(each_train_name){
      data.frame(m1 = each_name,m2 = each_train_name)
    }) %>%
      purrr::list_rbind()
  }) %>%
    purrr::list_rbind() %>%
    dplyr::filter(m1 != m2) %>%
    dplyr::filter(!(str_detect(m1,"^step") & str_detect(m2,"^step"))) %>%
    dplyr::filter(!(str_detect(m1,"^eNet|lasso|ridge") & str_detect(m2,"^eNet|lasso|ridge")))
}




