#' Training data using lasso-logistic
#'
#' @param input_df input_df
#' @param input_variables input_variables
#' @param input_y (0,1)
#' @param lambda lambda
#' @param seed seed
#' @param return_fit whether return fit
#'
#' @return model_coef, or list with fit
#' @export
#'
#' @examples
fun_train_logit <- function(input_df,input_variables,input_y,lambda="lambda.min",seed=1,return_fit=F){
  # 1)
  input_df %>%
    dplyr::select(dplyr::all_of(c(input_variables,input_y))) %>%
    na.omit()

  # 2) get x
  x = input_df %>%
    dplyr::select(dplyr::all_of(input_variables)) %>%
    as.matrix()
  y = input_df[[input_y]]

  set.seed(seed)
  fit = glmnet::cv.glmnet(x = as.matrix(x),y = y,family = "binomial")
  coef_data = glmnet::coef.glmnet(fit,s=lambda)
  variable_index = coef_data@i + 1
  variables_final = coef_data@Dimnames[[1]][variable_index]
  variables_value = coef_data@x

  coef_final =
    data.frame(row.names = variables_final,
               coef = variables_value)

  if(return_fit){
    list(fit = fit,
         coef = coef_final %>%
           dplyr::filter(rownames(.) %in% input_variables))
  }else{
    coef_final%>%
      dplyr::filter(rownames(.) %in% input_variables)
  }
}

