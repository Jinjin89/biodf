#' train data with randomForest
#'
#' @param input_df input_df
#' @param input_variables input_variables
#' @param input_y input_y, numeric for regression, factors for classification
#' @param seed seed
#' @param return_fit return_fit
#' @param top_feature The most important feature, default top 10
#' @param importance_method important methods, MeanDecreaseGini or MeanDecreaseAccuracy
#' @param randomForest_args other prams passing randomForst function
#'
#' @return train results of list or feature names
#' @export
#'
#' @examples
fun_train_logit_rf <-
  function (input_df, input_variables, input_y,seed = 1, return_fit = F,top_feature = 10,importance_method = "MeanDecreaseGini",randomForest_args = list(importance = TRUE)){
    input_df %>% dplyr::select(dplyr::all_of(c(input_variables,
                                               input_y))) %>% na.omit()
    x = input_df %>% dplyr::select(dplyr::all_of(input_variables)) %>%
      as.matrix()
    y = input_df[[input_y]]

    # 1). set random seed
    set.seed(seed)
    randomForest_args$x = x
    randomForest_args$y = y

    fit = do.call(randomForest::randomForest,randomForest_args)

    # 2) get importrance
    importrance <- fit$importance %>%
      as.data.frame() %>%
      dplyr::select(dplyr::all_of(importance_method)) %>%
      dplyr::arrange(dplyr::desc(!!as.name(importance_method))) %>%
      head(top_feature) %>%
      rownames()
    # 2).
    if (return_fit) {
      list(fit = fit,
           importrance = importrance)
    }else{
      importrance
    }
  }




#' logistic regression
#'
#' @param input_df input_df
#' @param input_variables input_variables
#' @param input_y input_y
#' @param fun_logit_multi_args fun_logit_multi_args
#'
#' @return data.frame
#' @export
#'
fun_logit_uni <- \(input_df, input_variables, input_y,
                   fun_logit_multi_args = list(return_pval = T)){
  # current_df
  input_df %<>% dplyr::select(dplyr::all_of(c(input_y,input_variables))) %>% na.omit
  fun_logit_multi_args$input_df <- input_df
  fun_logit_multi_args$input_y <- input_y

  purrr::map(input_variables,\(each_var){
    fun_logit_multi_args$input_variables = each_var
    do.call(fun_logit_multi,args = fun_logit_multi_args)
  }) %>%
    purrr::list_rbind()
}

#' logistic-regression: multiple
#'
#' @param input_df input_df
#' @param input_variables input_variables
#' @param input_y input_y
#' @param return_pval T
#' @param return_fit F
#' @param glm_args glm_args
#'
#' @return data.frame or list with fit and dataframe
#' @export
#'
fun_logit_multi <- \(input_df, input_variables, input_y,
                     return_pval=F,
                     return_fit=F,
                     glm_args = list()){
  input_df %<>% dplyr::select(dplyr::all_of(c(input_y,input_variables))) %>% na.omit
  # Construct the formula for logistic regression
  formula <- as.formula(paste(input_y, "~", paste(input_variables, collapse = "+")))

  # Fit logistic regression model
  glm_args$formula = formula
  glm_args$data = input_df
  glm_args$family = data = binomial

  #model <- stats::glm(formula, data = input_df, family = binomial)
  model <- do.call(stats::glm,glm_args)
  model_coef <- coef(model) %>%
    as.data.frame() %>%
    magrittr::set_colnames("coef") %>%
    dplyr::filter(rownames(.) != '(Intercept)')

  if(return_pval){
    # get confinit
    suppressMessages(
      {
        conf_int = stats::confint(model) %>%
          as.data.frame() %>%
          set_colnames(c("ci_low","ci_up"))
      }
    )

    coef(summary(model)) %>%
      as.data.frame() %>%
      dplyr::filter(rownames(.) != '(Intercept)') %>%
      dplyr::mutate(
        features = rownames(.),
        estimates = Estimate,
        pval = `Pr(>|z|)`,
        ci_low = conf_int[rownames(.),'ci_low'],
        ci_up = conf_int[rownames(.),'ci_up']
      ) %>%
      dplyr::select(features,estimates,ci_low,ci_up,pval) %>%
      return()
  }else if(return_fit){
    return(list(model_coef = model_coef,
                fit = model))
  }else{
    return(model_coef)
  }

}

#' glment train
#'
#' @param input_df input_df
#' @param input_variables input_variables
#' @param input_y input_y
#' @param lambda lambda.min or lambda.1se
#' @param alpha 0-1
#' @param cv.glment_args cv.glment_args
#' @param seed seed
#' @param return_fit return_fit
#'
#' @return data.frame or list with fit and data.frame
#' @export
#'

fun_train_logit <- function (input_df,
                             input_variables,
                             input_y,
                             lambda = "lambda.min",
                             alpha = 1,
                             cv.glment_args = list(),
                             seed = 1, return_fit = F) {
  input_df %<>% dplyr::select(dplyr::all_of(c(input_y,input_variables))) %>% na.omit

  x = input_df %>% dplyr::select(dplyr::all_of(input_variables)) %>%
    as.matrix()
  y = input_df[[input_y]]

  # fit step
  set.seed(seed)
  cv.glment_args$x = as.matrix(x)
  cv.glment_args$y = y
  cv.glment_args$family = binomial
  cv.glment_args$alpha = alpha
  fit = do.call(glmnet::cv.glmnet,cv.glment_args)

  # get coef
  coef_data = glmnet::coef.glmnet(fit, s = lambda)
  variable_index = coef_data@i + 1
  variables_final = coef_data@Dimnames[[1]][variable_index]
  variables_value = coef_data@x
  coef_final = data.frame(row.names = variables_final, coef = variables_value)
  if (return_fit) {
    list(fit = fit, coef = coef_final %>% dplyr::filter(rownames(.) %in%
                                                          input_variables))
  }
  else {
    coef_final %>% dplyr::filter(rownames(.) %in% input_variables)
  }
}



#' logistic-regression: step
#'
#' @param input_df input_df
#' @param input_variables input_variables
#' @param input_y input_y
#' @param direction direction
#' @param return_fit return_fit
#'
#' @return data.frame or list with fit and data.frame and
#' @export
#'
fun_train_logit_step <-  function(input_df,input_variables, input_y,direction,return_fit = F) {
  input_df %<>% dplyr::select(dplyr::all_of(c(input_y,input_variables))) %>% na.omit
  # Construct the formula for logistic regression
  formula <- as.formula(paste(input_y, "~", paste(input_variables, collapse = "+")))

  # Fit logistic regression model with stepwise selection
  set.seed(2)
  model <- stats::step(stats::glm(formula, data = input_df, family = binomial), direction = direction,trace=0)

  # Return the fitted model
  model_coef = coef(model) %>%
    as.data.frame() %>%
    magrittr::set_colnames("coef") %>%
    dplyr::filter(rownames(.) != "(Intercept)") %>%
    magrittr::set_rownames(stringr::str_remove_all(rownames(.),'`'))
  if(return_fit){
    list(model_coef =model_coef,
         fit = model)
  }else{
    return(model_coef)
  }
}


#' logistic-regression: gbm
#'
#' @param input_df input_df
#' @param input_variables input_variables
#' @param input_y input_y
#' @param seed seed
#' @param return_fit return_fit
#' @param distribution the distribution of y
#' @param n.trees how many trees, default 100
#' @param interaction.depth interaction.depth
#' @param cv.folds cv: defulat 10
#' @param gbm_args #' @param cv.folds
#'
#' @return vector or gbm object
#' @export
#'
fun_train_logit_gbm <- \(
  input_df,
  input_variables,
  input_y,
  seed = 1,
  return_fit=F,
  distribution = "bernoulli",
  n.trees = 100,
  interaction.depth = 1,
  cv.folds = 10,
  gbm_args = list()){

  library(gbm)
  # 1)
  input_df %<>% dplyr::select(dplyr::all_of(c(input_y,input_variables))) %>% na.omit

  formula <- as.formula(paste(input_y, "~", paste(input_variables, collapse = "+")))
  gbm_args$formula = formula
  gbm_args$data = input_df
  gbm_args$distribution = distribution
  gbm_args$n.trees = n.trees
  gbm_args$interaction.depth = interaction.depth
  gbm_args$cv.folds = cv.folds


  # Fit boosted logistic regression model
  set.seed(seed)
  model <- do.call(gbm::gbm, gbm_args)

  # get kepts genes
  n_trees_perf = gbm::gbm.perf(model)
  genes_kept <-
    gbm::relative.influence(model,n_trees_perf) %>%
    as.data.frame() %>%
    magrittr::set_colnames("imp") %>%
    dplyr::filter(imp > 0) %>%
    rownames()

  # Return the fitted model
  if(return_fit){
    list(genes_kept = genes_kept,
         fit = model)
  }else{
    return(genes_kept)
  }
}

#' svm-regression
#'
#' @param input_df input_df
#' @param input_variables input_variables
#' @param input_y input_y
#' @param seed seed
#' @param cost cost: default 1
#' @param scale scale: F
#' @param kernel default: linear, consider radial
#' @param return_fit T
#' @param svm_args smv_args
#'
#' @return data.frame or list with fit and data.frame
#' @export
#'

fun_train_logit_svm <- \(input_df,
                         input_variables,
                         input_y,
                         seed = 1,
                         cost = 1,
                         scale=F,
                         kernel = 'linear',
                         return_fit=F,
                         svm_args = list()){
  library(e1071)
  input_df %<>% dplyr::select(dplyr::all_of(c(input_y,input_variables))) %>% na.omit
  formula <- as.formula(paste(input_y, "~", paste(input_variables, collapse = "+")))
  final_args <- list()
  final_args$formula = formula
  final_args$data = input_df
  final_args$kernel = kernel
  final_args$cost = cost
  final_args$scale = scale
  for(each_name in names(svm_args)){
    final_args[[each_name]] <- svm_args[[each_name]]
  }
  # Fit SVM model
  model <- do.call(e1071::svm,final_args)

  # get the support vector index
  if(kernel == "linear"){
    model_coef = coef(model) %>%
      as.data.frame() %>%
      magrittr::set_colnames("coef") %>%
      dplyr::filter(rownames(.) != '(Intercept)')

  }else{
    model_coef = NULL
  }
  if(return_fit){
    list(model_coef = model_coef,
         fit = model) %>%
      return
  }else{
    return(model_coef)

  }
}


#' logistic feature selection step
#'
#' @param input_df input_df
#' @param input_variables input_variables
#' @param input_y input_y
#' @param outfile outfile
#' @param pval_cutoff pval_cutoff
#' @param methods_remove methods_remove
#' @param methods methods
#' @param fun_train_logit_args fun_train_logit_args
#' @param fun_train_logit_step_args fun_train_logit_step_args
#' @param fun_logit_multi_args fun_logit_multi_args
#' @param fun_logit_uni_args fun_logit_uni_args
#' @param fun_train_logit_gbm_args fun_train_logit_gbm_args
#' @param fun_train_logit_rf_args fun_train_logit_rf_args
#'
#' @return list with kept featurs
#' @export
#'
fun_feature_selection_logit <-function (
    input_df, input_variables, input_y,
    outfile,
    pval_cutoff = 0.05,
    methods_remove = NULL,
    methods = c("lasso", "eNet75", "eNet50","eNet25",
                "stepForward","stepBackward", "stepBoth",
                "multi_logit","uni_logit",
                "gbm", "rf"),
    fun_train_logit_args = list(),
    fun_train_logit_step_args = list(),
    fun_logit_multi_args = list(),
    fun_logit_uni_args = list(),
    fun_train_logit_gbm_args = list(),
    fun_train_logit_rf_args = list())
{
  if (!file.exists(outfile)) {
    message("results not found, running...")
    df_used <- input_df %>%
      dplyr::select(all_of(c(input_variables,input_y)))
    methods_select <- methods %>% dplyr::setdiff(methods_remove)
    message(
      "Current select feature selection methods: ",
      length(methods_select),"\n",
      paste0(methods_select, collapse = ","),
      "\nCurrent select feature selection features:",
      length(input_variables))
    # get args -> fun_train_logit_args
    fun_train_logit_args$input_df = input_df
    fun_train_logit_args$input_variables = input_variables
    fun_train_logit_args$input_y = input_y

    fun_train_logit_step_args$input_df = input_df
    fun_train_logit_step_args$input_variables = input_variables
    fun_train_logit_step_args$input_y = input_y

    fun_logit_multi_args$input_df = input_df
    fun_logit_multi_args$input_variables = input_variables
    fun_logit_multi_args$input_y = input_y
    fun_logit_multi_args$return_pval =T

    fun_logit_uni_args$input_df = input_df
    fun_logit_uni_args$input_variables = input_variables
    fun_logit_uni_args$input_y = input_y

    fun_train_logit_gbm_args$input_df = input_df
    fun_train_logit_gbm_args$input_variables = input_variables
    fun_train_logit_gbm_args$input_y = input_y

    fun_train_logit_rf_args$input_df = input_df
    fun_train_logit_rf_args$input_variables = input_variables
    fun_train_logit_rf_args$input_y = input_y

    tictoc::tic()
    results <- purrr::map(methods_select, function(each_method) {
      message("-------><-------")
      message(each_method, ": runing...")
      if (each_method == "lasso") {
        fun_train_logit_args$alpha = 1
        features_remained <- do.call(fun_train_logit, fun_train_logit_args) %>% rownames()
      }
      else if (each_method == "eNet75") {
        fun_train_logit_args$alpha = 0.75
        features_remained <- do.call(fun_train_logit, fun_train_logit_args) %>% rownames()
      }
      else if (each_method == "eNet50") {
        fun_train_logit_args$alpha = 0.5
        features_remained <- do.call(fun_train_logit, fun_train_logit_args) %>% rownames()
      }
      else if (each_method == "eNet25") {
        fun_train_logit_args$alpha = 0.25
        features_remained <- do.call(fun_train_logit, fun_train_logit_args) %>% rownames()
      }
      else if (each_method == "stepForward") {
        fun_train_logit_step_args$direction = "forward"
        features_remained <- do.call(fun_train_logit_step, fun_train_logit_step_args) %>% rownames()
      }
      else if (each_method == "stepBackward") {
        fun_train_logit_step_args$direction = "backward"
        features_remained <- do.call(fun_train_logit_step, fun_train_logit_step_args) %>% rownames()
      }
      else if (each_method == "stepBoth") {
        fun_train_logit_step_args$direction = "both"
        features_remained <- do.call(fun_train_logit_step, fun_train_logit_step_args) %>% rownames()
      }
      else if (each_method == "multi_logit") {
        message("cutoff: ", pval_cutoff)
        features_remained <- do.call(
          fun_logit_multi,
          fun_logit_multi_args) %>% dplyr::filter(
            pval <
              pval_cutoff) %>% rownames()
      }
      else if (each_method == "uni_logit") {
        message("cutoff: ", pval_cutoff)
        features_remained <- do.call(
          fun_logit_uni,
          fun_logit_uni_args) %>% dplyr::filter(
            pval <
              pval_cutoff) %>% rownames()
      }
      else if (each_method == "gbm") {
        features_remained <- do.call(fun_train_logit_gbm, fun_train_logit_gbm_args)
      }
      else if (each_method == "rf") {
        fun_train_logit_rf_args$input_df[[fun_train_logit_rf_args$input_y]] %<>%
          as.factor(.)
        features_remained <- do.call(fun_train_logit_rf, fun_train_logit_rf_args)
      }
      else {
        warning("invalida input: continue next:")
        warning("This is an invalid method input,check the input,or remove this methods")
        features_remained = "This is an invalid method input,check the input,or remove this methods"
      }
      message(each_method, "kept features: ", length(features_remained))
      return(features_remained)
    }) %>% set_names(methods_select)
    tictoc::toc()
    results %>% saveRDS(outfile)
  }
  message("Loading results")
  results = readRDS(outfile)
  for (each_method in names(results)) {
    message(each_method, ": ", length(results[[each_method]]))
  }
  return(results)
}






#' logit-train, with multiple methods
#'
#' @param input_df
#' @param input_y
#' @param design_df
#' @param input_features_list
#' @param outfile
#' @param feature_selection_method
#' @param train_method
#' @param fun_train_logit_args
#' @param fun_train_logit_step_args
#' @param fun_logit_multi_args
#' @param fun_train_logit_svm_args
#'
#' @return list of models
#' @export
#'

fun_train_multiple_logit <- function (input_df, input_y,design_df,
                                      input_features_list, outfile,
                                      feature_selection_method = "m1",
                                      train_method = "m2",
                                      fun_train_logit_args = list(),
                                      fun_train_logit_step_args = list(),
                                      fun_logit_multi_args = list(),
                                      fun_train_logit_svm_args = list()) {
  if (!file.exists(outfile)) {
    message("model not found, running...")
    results = list()
    input_features_list = purrr::map(input_features_list,
                                     function(each_var) {
                                       each_var %>% stringr::str_remove_all("`")
                                     })
    all_var_check = input_features_list %>% unlist() %>%
      unname() %>% unique()
    df_check <- input_df %>% dplyr::select(dplyr::all_of(c(input_y,all_var_check)))

    fun_train_logit_args$input_df = input_df
    fun_train_logit_args$input_y = input_y

    fun_train_logit_step_args$input_df = input_df
    fun_train_logit_step_args$input_y = input_y

    fun_logit_multi_args$input_df = input_df
    fun_logit_multi_args$input_y = input_y

    fun_train_logit_svm_args$input_df = input_df
    fun_train_logit_svm_args$input_y = input_y

    for (i in seq_along(design_df[[1]])) {
      m1 = design_df[[feature_selection_method]][i]
      each_method = design_df[[train_method]][i]
      model_name = paste0(m1, "&", each_method)
      message("-------------><-------------")
      message("model_name: ", model_name)
      input_variables = input_features_list[[m1]]
      fun_train_logit_args$input_variables = input_variables
      fun_train_logit_step_args$input_variables = input_variables
      fun_logit_multi_args$input_variables = input_variables
      fun_train_logit_svm_args$input_variables = input_variables

      if (each_method == "lasso") {
        fun_train_logit_args$alpha = 1
        model_df <- do.call(fun_train_logit, fun_train_logit_args)
      }
      else if (each_method == "eNet75") {
        fun_train_logit_args$alpha = 0.75
        model_df <- do.call(fun_train_logit, fun_train_logit_args)
      }
      else if (each_method == "eNet50") {
        fun_train_logit_args$alpha = 0.5
        model_df <- do.call(fun_train_logit, fun_train_logit_args)
      }
      else if (each_method == "eNet25") {
        fun_train_logit_args$alpha = 0.25
        model_df <- do.call(fun_train_logit, fun_train_logit_args)
      }
      else if (each_method == "ridge") {
        fun_train_logit_args$alpha = 0
        model_df <- do.call(fun_train_logit, fun_train_logit_args)
      }
      else if (each_method == "stepForward") {
        fun_train_logit_step_args$direction = "forward"
        model_df <- do.call(fun_train_logit_step, fun_train_logit_step_args)
      }
      else if (each_method == "stepBackward") {
        fun_train_logit_step_args$direction = "backward"
        model_df <- do.call(fun_train_logit_step, fun_train_logit_step_args)
      }
      else if (each_method == "stepBoth") {
        fun_train_logit_step_args$direction = "both"
        model_df <- do.call(fun_train_logit_step, fun_train_logit_step_args)
      }
      else if (each_method == "multi_logit") {
        model_df <- do.call(fun_logit_multi, fun_logit_multi_args)
      }
      else if (each_method == "svm") {
        model_df <- do.call(fun_train_logit_svm,fun_train_logit_svm_args)
      }
      else {
        warning("invalid input: continue next:")
        warning("This is an invalid method input,check the input,or remove this methods")
        model_df = "This is an invalid method input,check the input,or remove this methods"
      }
      message("model variables count", ": ", nrow(model_df))
      message("model variables: ", paste0(rownames(model_df),
                                          collapse = ","))
      results[[model_name]] = model_df
    }
    results %>% saveRDS(outfile)
  }
  results = readRDS(outfile)

  for (model_name in names(results)) {
    model_df = results[[model_name]]
    message("model_name", ": ", model_name)
    message("model variables count", ": ", nrow(model_df))
    message("model variables: ",
            paste0(rownames(model_df),
                   collapse = ","))
  }
  results
}
