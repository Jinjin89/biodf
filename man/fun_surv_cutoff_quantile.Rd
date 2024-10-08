% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/fun_surv.R
\name{fun_surv_cutoff_quantile}
\alias{fun_surv_cutoff_quantile}
\title{cutoff with quantile points, e.g. if you supplied input_pct as c(0.25,0.5,0.75), you will get four groups}
\usage{
fun_surv_cutoff_quantile(
  input_df,
  input_variables = "Risk_score",
  output_variables_name = "Group4",
  input_pct = c(0.25, 0.5, 0.75),
  input_labels = c("Low", "Intermediate_low", "Intermediate_high", "High"),
  factor_level = T
)
}
\arguments{
\item{input_df}{input_df}

\item{input_variables}{input_variables}

\item{output_variables_name}{output_variables_name}

\item{input_pct}{perentage should be in (0,1)}

\item{input_labels}{the lables length should be length(input_pct)+1}

\item{factor_level}{T or F, whether to refactor as input_labels}
}
\value{
todo
}
\description{
cutoff with quantile points, e.g. if you supplied input_pct as c(0.25,0.5,0.75), you will get four groups
}
