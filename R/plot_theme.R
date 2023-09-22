
#' Default theme
#'
#' @param ... passing to theme
#'
#' @return theme object
#' @export


tn = function(...){
  ggplot2::theme_minimal()+
    ggplot2::theme(
      panel.grid = ggplot2::element_blank(),
      axis.line = ggplot2::element_line(),
      axis.ticks = ggplot2::element_line())+
    ggplot2::theme(...)
}

#' theme for barplot
#'
#' @param ... passing to theme
#' @param show_axis_title  whether plot axis title
#' @param legend_text  legend text size
#'
#' @return theme obj
#' @export
#'
tn_bar = function(show_axis_title=F,legend_text = 6,...){
  ggplot2::theme_minimal()+
    ggplot2::theme(
      panel.grid = ggplot2::element_blank(),
      panel.background = ggplot2::element_rect(fill =NA,color = "black"),
      axis.ticks = ggplot2::element_line(color = "black"),
      axis.title = element_blank(),
      legend.title = ggplot2::element_blank(),
      legend.text = ggplot2::element_text(size = legend_text),
      legend.key.size = ggplot2::unit(4,"mm"),
      legend.position = "top"
    )+
    ggplot2::theme(...) -> tn

  if(show_axis_title){
    tn =
      tn + ggplot2::theme(
        axis.title = element_text()
      )
  }
  tn
}


#' No legend
#'
#' @param ... pass to ggplot theme function
#'
#' @return
#' @export
#'
tn_no_legend <- function(...){
  ggplot2::theme(legend.position = "none",...)
}

