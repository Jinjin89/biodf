% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/fun_surv.R
\name{fun_plot_surv_cutoff}
\alias{fun_plot_surv_cutoff}
\title{Plot unicox pvalues and HR with differnet cutoff}
\usage{
fun_plot_surv_cutoff(
  input_df,
  input_variables = "Risk_score",
  input_pct = seq(0.2, 0.8, 0.01),
  axis_fold = 2,
  return_cutoff = F,
  x_label = "Quantile(\%)",
  y_label = "P-value",
  sec_y_axis_label = "Hazard ratio",
  grid_color = "gray90",
  vline_lty = "twodash",
  vline_color = "#FF9B82",
  line_color = c("#8BE8E5", "#9400FF")
)
}
\arguments{
\item{input_df}{input_df}

\item{input_variables}{input_variables}

\item{input_pct}{input_pct,deafulat from 0.2 to 0.8 by 0.01}

\item{axis_fold}{axis_fold is 2}

\item{x_label}{x_label}

\item{y_label}{y_label}

\item{sec_y_axis_label}{sec_y_axis_label}

\item{grid_color}{grid_color}

\item{vline_lty}{vline_lty}

\item{vline_color}{vline_color}

\item{line_color}{line_color}
}
\value{
todo
}
\description{
Plot unicox pvalues and HR with differnet cutoff
}
