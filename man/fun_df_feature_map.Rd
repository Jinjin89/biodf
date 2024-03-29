% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/df.R
\name{fun_df_feature_map}
\alias{fun_df_feature_map}
\title{relabel features with new label in element-wise way}
\usage{
fun_df_feature_map(
  input_df,
  input_variables,
  from,
  to,
  new_name = NULL,
  not_found_replace = NA_character_
)
}
\arguments{
\item{input_df}{data.frame}

\item{input_variables}{the variable name in data.frame}

\item{from}{the old element label}

\item{to}{the new element label}

\item{new_name}{new_name, if not supplied, change inplace}

\item{not_found_replace}{the label}
}
\value{
relabeled dataframe for supplied variable
}
\description{
relabel features with new label in element-wise way
}
\examples{
fun_df_feature_map(biodata::tcga_clin_xena,
  input_variables = "gender",
  from = c("MALE","FEMALE"),
  to = c("M","F"),
  new_name = "sex")
}
