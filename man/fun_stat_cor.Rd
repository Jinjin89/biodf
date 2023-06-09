% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/stat_compare.R
\name{fun_stat_cor}
\alias{fun_stat_cor}
\title{Calculate the correlation between the input x and y one by one}
\usage{
fun_stat_cor(input_df, input_xs, input_ys, method = "sp")
}
\arguments{
\item{input_df}{the input data.frame}

\item{input_xs}{input x ,should be in input_df}

\item{input_ys}{input y ,should be in input_df}

\item{method}{spearman or person}
}
\value{
data.frame of the correlation results
}
\description{
Calculate the correlation between the input x and y one by one
}
