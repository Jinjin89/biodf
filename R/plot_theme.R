#' Default theme
#'
#' @param ... passing to theme
#'
#' @return theme object
#' @export


tn = function(...){
  ggplot2::theme_minimal()+
    ggplot2::theme(
      text = ggplot2::element_text(),
      panel.grid = ggplot2::element_blank(),
      #plot.background = ggplot2::element_blank(),
      plot.background = ggplot2::element_rect(fill = NA,colour = NA),
      legend.background = ggplot2::element_blank(),
      panel.background = ggplot2::element_rect(),
      axis.line = ggplot2::element_blank(),
      axis.ticks = ggplot2::element_line(),
      ...)
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
tn_bar = function(show_axis_title=F,show_legend_title=F,legend_text = 6,...){
  ggplot2::theme_minimal()+
    ggplot2::theme(
      text = ggplot2::element_text(),
      panel.grid = ggplot2::element_blank(),
      panel.background = ggplot2::element_rect(fill =NA,color = "black"),
      axis.ticks = ggplot2::element_line(color = "black"),
      axis.title = element_blank(),
      legend.title = ggplot2::element_blank(),
      legend.text = ggplot2::element_text(size = legend_text),
      legend.key.size = ggplot2::unit(4,"mm"),
      legend.position = "top",
      ...
    ) -> tn

  if(show_axis_title){
    tn =
      tn + ggplot2::theme(
        axis.title = element_text()
      )
  }
  if(show_legend_title){
    tn = tn+ggplot2::theme(
      legend.title = ggplot2::element_text()
    )
  }
  tn
}


#' No legend
#'
#' @param ... pass to ggplot theme function
#'
#' @return a ggplot theme
#' @export
#'
tn_no_legend <- function(...){
  ggplot2::theme(legend.position = "none",...)
}

#' theme empty
#'
#' @param ...
#'
#' @return ggplot-theme
#' @export
#'
tn_empty <- function(...){
  theme_minimal()+
    theme(axis.title = element_blank())+
    theme(panel.grid = element_blank())+
    theme(...)
}


#' italic the text x
#'
#' @param ...  passing to element_text
#'
#' @return ggplot2 theme
#' @export
#'
tn_italic_text_x <- function(...){
  ggplot2::theme(
    axis.text.x = ggplot2::element_text(face = 'italic',...)
  )
}


#' italic the text y
#'
#' @param ...  passing to element_text
#'
#' @return ggplot2 theme
#'
#'
tn_italic_text_y <- function(...){
  ggplot2::theme(
    axis.text.y= ggplot2::element_text(face = 'italic',...)
  )
  }

#' italic the title x
#'
#' @param ...  passing to element_text
#'
#' @return ggplot2 theme
#'
tn_italic_title_x <- function(...){
  ggplot2::theme(
    axis.title.x = ggplot2::element_text(face = 'italic',...)
  )
}

#' italic the title y
#'
#' @param ...  passing to element_text
#'
#' @return ggplot2 theme
#'
tn_italic_title_y <- function(...){
  ggplot2::theme(
    axis.title.y= ggplot2::element_text(face = 'italic',...)
  )}
