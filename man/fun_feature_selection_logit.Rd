% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/logistic.R
\name{fun_feature_selection_logit}
\alias{fun_feature_selection_logit}
\title{logistic feature selection step}
\usage{
fun_feature_selection_logit(
  input_df,
  input_variables,
  input_y,
  outfile,
  pval_cutoff = 0.05,
  methods_remove = NULL,
  methods = c("lasso", "eNet75", "eNet50", "eNet25", "stepForward", "stepBackward",
    "stepBoth", "multi_logit", "uni_logit", "gbm", "rf"),
  fun_train_logit_args = list(),
  fun_train_logit_step_args = list(),
  fun_logit_multi_args = list(),
  fun_logit_uni_args = list(),
  fun_train_logit_gbm_args = list(),
  fun_train_logit_rf_args = list()
)
}
\arguments{
\item{input_df}{input_df}

\item{input_variables}{input_variables}

\item{input_y}{input_y}

\item{outfile}{outfile}

\item{pval_cutoff}{pval_cutoff}

\item{methods_remove}{methods_remove}

\item{methods}{methods}

\item{fun_train_logit_args}{fun_train_logit_args}

\item{fun_train_logit_step_args}{fun_train_logit_step_args}

\item{fun_logit_multi_args}{fun_logit_multi_args}

\item{fun_logit_uni_args}{fun_logit_uni_args}

\item{fun_train_logit_gbm_args}{fun_train_logit_gbm_args}

\item{fun_train_logit_rf_args}{fun_train_logit_rf_args}
}
\value{
list with kept featurs
}
\description{
logistic feature selection step
}
