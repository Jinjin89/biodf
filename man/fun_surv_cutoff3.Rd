% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/fun_surv.R
\name{fun_surv_cutoff3}
\alias{fun_surv_cutoff3}
\title{group var into three groups}
\usage{
fun_surv_cutoff3(
  input_df,
  input_variables = "Risk_score",
  output_variables_name = "Group3",
  input_pct = c(0.25, 0.75),
  input_labels = c("Low", "Medium", "High"),
  factor_level = T
)
}
\value{
df
}
\description{
group var into three groups
}
