#' nomogram predict function
#'
#' @param input_df input data frame with each variables
#' @param input_fit the trained model
#' @param input_variables the variables used in fit
#' @param input_y the response data
#'
#' @return a predictive dataframe with yhat(which is the predictive results)
#' @export
#'

fun_nomogram_predict <- function(input_df,input_fit,input_variables,input_y){

  # 1) initialize
  input_args_each = list()
  input_args_each$object = input_fit
  input_args_each$fun = plogis

  print("Precition step")
  print("input_shape:")
  print(dim(input_df))
  input_df = na.omit(input_df)
  print("remove na shape:")
  print(dim(input_df))


  # 2) loop to predict each
  purrr::map_df(seq_along(input_df[[1]]),function(i){
    # 1) get args
    for(each_var in input_variables){
      input_args_each[[each_var]] = input_df[i,each_var]
    }
    out = do.call("Predict",input_args_each)
    as.data.frame(out)
  }) ->df_out

  # 3) merge label into preeicted results
  df_out$y = input_df[[input_y]]


  # 2)
  df_out
}

#' nomogram training for logistics regression
#'
#' @param input_df the training data
#' @param input_variables the training varaibles
#' @param input_y the response data,one-hot encoding is required
#' @param val_list the predictive results,after training, predict the yhat for each validation datasets
#' @param xfrac nomogram plot function
#' @param cex.axis nomogram plot function
#' @param cex.var nomogtam plot funciton
#'
#' @return todo
#' @export
#'

fun_nomogram_logit <- function(input_df,input_variables,input_y,
                               val_list = NULL,
                               funlabel = "Probability",
                               xfrac=0.2,
                               cex.axis=.5,
                               cex.var=0.8){

  require(rms)
  print("input_shape:")
  print(dim(input_df))
  input_df %<>%
    select(all_of(c(input_variables,input_y))) %>%
    na.omit()
  print("remove na shape:")
  print(dim(input_df))

  fun_utils_assgin_global("dd_for_the_tmp",rms::datadist(input_df))
  options(datadist='dd_for_the_tmp')

  # 1) fit a logit
  f <- rms::lrm(
    as.formula(paste0(input_y,"~.")),
    data = input_df)
  # 2)
  nom <- rms::nomogram(f,fun=plogis, funlabel=funlabel)
  plot(nom, xfrac=xfrac,cex.axis=cex.axis,cex.var=cex.var)

  # 3) remove dd
  rm(dd_for_the_tmp,envir = .GlobalEnv)

  # 4) prediction
  df_out <- fun_nomogram_predict(
    input_df = input_df,
    input_fit = f,
    input_variables = input_variables,
    input_y = input_y
  )
  if(is.null(val_list)){
    return(df_out)
  }else{
    return_list = vector("list",length(val_list) + 1)
    return_list[[1]] = df_out
    for(i in seq_along(val_list)){
      return_list[[i+1]] =
        fun_nomogram_predict(
          input_df = val_list[[i]],
          input_fit = f,
          input_variables = input_variables,
          input_y = input_y
        )
    }
    return(return_list)
  }

}
