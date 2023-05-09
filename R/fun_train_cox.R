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
  coef_final = as.data.frame(coef(step_cox))
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

  coef_final = as.data.frame(coef(cb_fit))
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
