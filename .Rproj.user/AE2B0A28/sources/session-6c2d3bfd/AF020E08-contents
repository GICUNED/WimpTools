
# WimpGrid Biplot ------------------------------------------------

#' Weigthed Implication Grid Biplot
#'
#' @description This function generates a biplot visualization from the
#'              semi-structured interview data collected via the WimpGrid
#'              methodology. A biplot is a graphical representation that
#'              combines both hipothetical selves and personal constructs of a
#'              matrix into a single plot, making it easier to identify
#'              patterns, relationships, and clusters within the data.
#'
#' @param wimp Subject's WimpGrid. It must be a "wimp" S3 object
#'        imported by the \code{\link{importwimp}} function.
#' @param text.size Size of the text labels. Default is 1.
#'
#' @return A interactive biplot made with plotly
#'
#' @importFrom DescTools CartToPol PolToCart
#' @import useful
#' @importFrom stats cor prcomp runif
#' @export
#'
#' @examples
#'
#' wimp_biplot(example.wimp)
#'

wimp_biplot <- function(wimp, text.size = 1){

  wimp <- .align.wimp(wimp, exclude.dilemmatics = FALSE)

  lpoles <- wimp$constructs$left.poles
  rpoles <- wimp$constructs$right.poles
  dim <- length(lpoles)

  matrix <- .hypo.matrix(wimp)

  pca <- prcomp(matrix, rank. = 2)

  s.pc1 <- round(summary(pca)[6]$importance[2,1] * 100, digits = 2)
  s.pc2 <- round(summary(pca)[6]$importance[2,2] * 100, digits = 2)

  pc1 <-pca$rotation[,1]
  pc2 <-pca$rotation[,2]

  pc.r <- DescTools::CartToPol(pc1,pc2)$r
  pc.t <- DescTools::CartToPol(pc1,pc2)$theta

  names <- colnames(matrix)

  v1 <- c(pca$x[,1], -pca$x[,1])
  v2 <- c(pca$x[,2], -pca$x[,2])

  v.r <- max(abs(pc.r)) + 0.05 * max(abs(pc.r)) + max(abs(pc.r)) * runif(dim,-0.35,0.1)
  v.t <-  DescTools::CartToPol(v1,v2)$theta

  v1 <- DescTools::PolToCart(v.r,v.t)$x
  v2 <- DescTools::PolToCart(v.r,v.t)$y

  vnames <- c(rpoles,lpoles)

  df <- data.frame(names,pc1,pc2)
  dfv <- data.frame(vnames,v1,v2)

  range <- max(abs(dfv[2:3]) + 0.1 * abs(dfv[2:3]))

  fig <- plot_ly(type = "scatter") %>%
    add_annotations(
      data = df[-c(1,dim+2),],
      x = ~pc1,
      y = ~pc2,
      text = ~names,
      hoverinfo = 'text',
      font = ~list(size = 15 * text.size),
      showarrow = FALSE,
      xanchor = 'center',
      yanchor = 'center'
    ) %>%
    add_annotations(
      data = df[1,],
      x = ~pc1,
      y = ~pc2,
      text = ~names,
      hoverinfo = 'text',
      font = ~list(size = 15 * text.size, color="darkblue"),
      showarrow = FALSE,
      xanchor = 'center',
      yanchor = 'center'
    ) %>%
    add_annotations(
      data = df[dim+2,],
      x = ~pc1,
      y = ~pc2,
      text = ~names,
      hoverinfo = 'text',
      font = ~list(size = 15 * text.size, color="darkgreen"),
      showarrow = FALSE,
      xanchor = 'center',
      yanchor = 'center'
    ) %>%
    add_segments(
      data = dfv,
      x = 0, xend = ~v1,
      y = 0, yend = ~v2,
      line = list(color = '#FA9D13', dash = "dot", width = 0.75),
      hoverinfo = 'none',
      inherit = FALSE,
      showlegend = FALSE
    ) %>%
    add_annotations(
      data = dfv,
      x = ~v1,
      y = ~v2,
      text = ~vnames,
      hoverinfo = 'none',
      font = list(size = 12 * text.size, color = "#FA9D13"),
      showarrow = FALSE,
      xanchor = 'center',
      yanchor = 'bottom'
    ) %>%
    layout(
      xaxis = list(title = paste("<B>PC1</B> [", s.pc1, "%]", sep = ""),
                   range = c(-range, range),
                   gridcolor = "white",
                   gridwidth = 3,
                   zeroline = TRUE,
                   zerolinecolor = "black",
                   zerolinewidth = 2,
                   showline = FALSE),
      yaxis = list(title = paste("<B>PC2</B> [", s.pc2, "%]", sep = ""),
                   range = c(-range, range),
                   gridcolor = "white",
                   gridwidth = 3,
                   zeroline = TRUE,
                   zerolinecolor = "black",
                   zerolinewidth = 2,
                   showline = FALSE),
      showlegend = FALSE,
      plot_bgcolor = "#FFFAED"
    )

  return(fig)

}
