# Monitoring Adjustment Radar Plot ---------------------------------------------------

#' Monitoring Self Adjustment -- monitoring_adj()
#'
#' @description This function generates an interactive radar plot to visualize
#' and compare the progress of the "self" and "ideal" across different
#' time points (e.g., pre- and post-intervention).
#'
#' @param wimp.t0 First Subject's WimpGrid (e.g., pre-intervention). It must be a "wimp" S3 object.
#' @param wimp.t1 Second Subject's WimpGrid (e.g., post-intervention). It must be a "wimp" S3 object.
#' @param legend If TRUE, displays legend of plot. Default is TRUE.
#'
#' @return Interactive radar plot using plotly
#'
#' @import plotly
#' @export
#'
#' @examples
#'
#'  monitoring_adj(example.wimp,example.wimp)

monitoring_adj <- function(wimp.t0, wimp.t1, legend = TRUE){

  wimp.t0 <- .align.wimp(wimp.t0, exclude.dilemmatics = FALSE)
  wimp.t1 <- .align.wimp(wimp.t1, exclude.dilemmatics = FALSE)

  merge <- .merge.wimp(wimp.t0,wimp.t1)
  if(.compatibility.merge.wimp(wimp.t0,wimp.t1) == "Incompatibility"){stop("WimpGrids have no constructs in common. Monitoring not possible.")}

  self.t0 <- wimp.t0$self$standarized[merge$index1]
  self.t0 <- c(self.t0,self.t0[1])

  ideal.t0 <- wimp.t0$ideal$standarized[merge$index1]
  ideal.t0 <- c(ideal.t0,ideal.t0[1])

  r.poles.t0 <- wimp.t0$constructs$right.poles[merge$index1]
  l.poles.t0 <- wimp.t0$constructs$left.poles[merge$index1]
  poles.t0 <- paste(r.poles.t0," (",l.poles.t0,")", sep="")
  poles.t0 <- c(poles.t0,poles.t0[1])

  construct.t0 <- wimp.t0$constructs$constructs[merge$index1]
  construct.t0 <- c(construct.t0,construct.t0[1])

  colors.t0 <- .construct.colors(wimp.t0, mode = "red/green")[merge$index1,1]
  colors.t0 <- c(colors.t0,colors.t0[1])


  self.t1 <- wimp.t1$self$standarized[merge$index2]
  self.t1 <- c(self.t1,self.t1[1])

  ideal.t1 <- wimp.t1$ideal$standarized[merge$index2]
  ideal.t1 <- c(ideal.t1,ideal.t1[1])

  r.poles.t1 <- wimp.t1$constructs$right.poles[merge$index2]
  l.poles.t1 <- wimp.t1$constructs$left.poles[merge$index2]
  poles.t1 <- paste(r.poles.t1," (",l.poles.t1,")", sep="")
  poles.t1 <- c(poles.t1,poles.t1[1])

  construct.t1 <- wimp.t1$constructs$constructs[merge$index2]
  construct.t1 <- c(construct.t1,construct.t1[1])

  colors.t1 <- .construct.colors(wimp.t1, mode = "red/green")[merge$index2,1]
  colors.t1 <- c(colors.t1,colors.t1[1])

  fig <- plot_ly(
    type = 'scatterpolar',
    fill = 'toself'
  )
  fig <- fig %>%
    add_trace(
      mode = "lines",
      r = 0,
      theta = poles.t0,
      fill = "none",
      line = list(color = "#444444", width = 1.5, shape = 'spline', smoothing = 1),
      name = 'Pole Threshold',
      hoverinfo = 'none',
      showlegend = FALSE
    )
  fig <- fig %>%
    add_trace(
      mode = "lines",
      r = ideal.t0,
      theta = poles.t0,
      fill = "none",
      line = list(color = "darkgreen", width = 2.5, dash ="dot"),
      name = 'Baseline Ideal',
      hoverinfo = 'none'
    )
  fig <- fig %>%
    add_trace(
      mode = "lines",
      r = ideal.t1,
      theta = poles.t1,
      fill = "none",
      line = list(color = "darkgreen", width = 3),
      name = 'Actual Ideal',
      hoverinfo = 'none'
    )
  fig <- fig %>%
    add_trace(
      r = self.t0,
      theta = poles.t0,
      name = "Baseline",
      marker = list(color = colors.t0, size = 7, line = list(color = '#FA9D13', width = 1.5)),
      fillcolor = 'rgba(255, 217, 125, 0.5)',
      line = list(width = 1, color = "#FA9D13"),
      text = ~paste('<B>',construct.t0,'</B>', '\nSelf:', round(self.t0, 2), '\nIdeal:', round(ideal.t0,2)),
      hoverinfo = 'text',
      hoverlabel=list(bgcolor = colors)
    )

  fig <- fig %>%
    add_trace(
      r = self.t1,
      theta = poles.t1,
      name = "Actual",
      marker = list(color = colors.t1, size = 7, line = list(color = '#AB81E5', width = 1.5)),
      fillcolor = 'rgba(213, 192, 242, 0.5)',
      line = list(width = 1, color = "#AB81E5"),
      text = ~paste('<B>',construct.t1,'</B>', '\nSelf:', round(self.t1, 2), '\nIdeal:', round(ideal.t1,2)),
      hoverinfo = 'text',
      hoverlabel=list(bgcolor = colors)
    )
  fig <- fig %>%
    layout(
      showlegend = legend,
      polar = list(
        radialaxis = list(
          visible = T,
          range = c(-1,1)
        )
      )
    )

  fig

}

# Monitoring Heatmap -----------------------------------------------------------

#' Monitoring SSI Adjustment -- monitoring_ssi()
#'
#' @description This function generates an interactive heatmap to visualize
#' and compare the progress of the SSI index across different
#' time points (e.g., pre- and post-intervention).
#'
#' @param wimp.t0 First Subject's WimpGrid (e.g., pre-intervention). It must be a "wimp" S3 object.
#' @param wimp.t1 Second Subject's WimpGrid (e.g., post-intervention). It must be a "wimp" S3 object.
#'
#' @return A two heatmap made with plotly.
#'
#' @import plotly
#' @export
#'
#' @examples
#'
#' monitoring_ssi(example.wimp, example.wimp)
#'


monitoring_ssi <- function(wimp.t0, wimp.t1){

  create_heatmap <- function(wimp, show_y_axis_title = TRUE, show_legend = FALSE, hide_y_ticks = FALSE) {
    x <- wimp$self$standarized
    y <- wimp$ideal$standarized

    alpha.values <- seq(0, 1, by = 0.01)
    beta.values <- seq(0, 1, by = 0.01)

    sim_matrix <- outer(alpha.values, beta.values, Vectorize(function(alpha, beta) {
      .sim_index(x, y, alpha = alpha, beta = beta)
    }))

    plot_ly(
      x = alpha.values,
      y = beta.values,
      z = sim_matrix,
      type = "heatmap",
      colorscale = list(c(0, "#F52722"), c(0.5, "yellow"), c(1, "#A5D610")),
      zmin = 0,
      zmax = 1,
      hovertemplate = '<b>Alpha:</b> %{x}<br><b>Beta:</b> %{y}<br><b>Adjustment:</b> %{z}<extra></extra>',
      showscale = show_legend
    ) %>%
      layout(
        xaxis = list(title = "Attention to Self Discrepances"),
        yaxis = list(
          title = if (show_y_axis_title) "Attention to the Desired Change" else NULL,
          showticklabels = !hide_y_ticks,
          ticks = if (hide_y_ticks) "" else "outside"
        ),
        shapes = list(
          list(
            type = "rect",
            x0 = -0.005,
            x1 = 1.005,
            y0 = -0.005,
            y1 = 1.005,
            line = list(color = "black", width = 2)),
          list(
            type = "line",
            x0 = 0,
            y0 = 1,
            x1 = 1,
            y1 = 0,
            line = list(
              color = "black",
              width = 1,
              dash = "dot"
            )
          )
        ),
        annotations = list(
          list(
            x = 0.5,
            y = 0.5,
            xref = "x",
            yref = "y",
            text = "+",
            showarrow = FALSE,
            font = list(color = "black", size = 20)
          )
        )
      ) %>%
      style(
        hoverlabel = list(
          bgcolor = 'rgba(255, 255, 255, 0.8)',
          bordercolor = 'black',
          font = list(size = 12)
        )
      )
  }

  heatmap1 <- create_heatmap(wimp.t0, show_y_axis_title = TRUE, show_legend = TRUE) %>% layout(title = "Pre-Intervention")
  heatmap2 <- create_heatmap(wimp.t1, show_y_axis_title = FALSE, show_legend = FALSE, hide_y_ticks = TRUE) %>% layout(title = "Post-Intervention")

  subplot(heatmap1, heatmap2, nrows = 1, titleX = TRUE, titleY = TRUE, margin  = 0.005) %>%
    layout(
      title = "",
      showlegend = FALSE
    )
}

# Monitoring PH Index -------------------------------------------------------------------------------
#' Monitoring PH index -- monitoring_ph()
#'
#' @description This function generates a graphical comparison of constructs from two
#' different Weigthed Implications Grids plotted in a PH space.
#
#' @param wimp.t0 Data object containing constructs and their respective P and H coordinates for the first grid.
#' @param wimp.t1 Data object containing constructs and their respective P and H coordinates for the second grid
#' @param show.centroid Logical; if TRUE, displays the centroid of construct P-H coordinates for both grids on the graph.
#' @param text.size Size of the text labels. Default is 1.
#' @param ... additional arguments are passed from \code{\link{ph_index}} function.
#'
#' @return A Plotly object representing the comparative graph of constructs across both grids
#'
#' @export
#' @examples
#'
#' monitoring_ph(example.wimp,example.wimp)

monitoring_ph <- function(wimp.t0, wimp.t1, show.centroid = TRUE, text.size = 1
                          ,...) {

  wimp.t0 <- .align.wimp(wimp.t0, exclude.dilemmatics = FALSE)
  wimp.t1 <- .align.wimp(wimp.t1, exclude.dilemmatics = FALSE)

  merge <- .merge.wimp(wimp.t0,wimp.t1)
  if(.compatibility.merge.wimp(wimp.t0,wimp.t1) == "Incompatibility"){stop("WimpGrids have no constructs in common. Monitoring not possible.")}

  # Calculate the Mahalanobis distance matrix for wimp.t1
  phm.matII <- ph_index(wimp = wimp.t1, ...)[merge$index2,]
  phm.matII.df <- as.data.frame(phm.matII)
  phm.matII.df$constructo <- rownames(phm.matII)
  phm.matII.df$self.constr <- wimp.t1$constructs$self.poles[merge$index2]

  # Calculate the Mahalanobis distance matrix for wimp.t0
  phm.matI <- ph_index(wimp = wimp.t0, ...)[merge$index1,]
  phm.matI.df <- as.data.frame(phm.matI)
  phm.matI.df$constructo <- rownames(phm.matI)
  phm.matI.df$self.constr <- wimp.t0$constructs$self.poles[merge$index1]

  # Define the boundaries of the regions according to wimp.graphII
  limit <- max(abs(phm.matII.df$p), abs(phm.matII.df$h)) * 1.1

  # Opacity for wimp.graphI constructs
  wgI.opacity <- 0.5
  wgI.color <- 'grey' #gray(0.83)

  # Configuring regions using phm.matII
  shapes <- list(
    list(type = "path", path = paste("M 0,0 L", limit, ",", limit, " L0,", limit, " Z"),
         fillcolor = "#FFD97D", opacity = 0.2, line = list(color = "#FA9D13")),
    list(type = "path", path = paste("M 0,0 L", limit, ",", -limit, " L0,", -limit, " Z"),
         fillcolor = "#FFD97D", opacity = 0.2, line = list(color = "#FA9D13")),
    list(type = "line", x0 = 0, y0 = 0, x1 = limit, y1 = limit,
         xref = "x", yref = "y", line = list(color = "#FFD97D", width = 1, dash = "dash")),
    list(type = "line", x0 = 0, y0 = 0, x1 = limit, y1 = -limit,
         xref = "x", yref = "y", line = list(color = "#FFD97D", width = 1, dash = "dash"))
  )

  # Initialising the Plotly chart
  p <- plot_ly()

  # Configuring layout
  p <- p %>%
    layout(title = '',
           xaxis = list(title = 'Presence'),
           yaxis = list(title = 'Hierarchy'),
           plot_bgcolor = "white",
           font = list(family = "Arial"),
           showlegend = FALSE,
           shapes = shapes)

  # Add wimp.graphI dots (with lower opacity)
  colorsI <- .construct.colors(wimp = wimp.t0, mode = "red/green")
  phm.matI.df$color <- colorsI[merge$index1,"color"]
  p <- p %>%
    add_markers(data = phm.matI.df, x = ~p, y = ~h,
                marker = list(color = ~color, size = 4, opacity = wgI.opacity,
                              line = list(color = 'black', width = 1, dash = "dot")),
                text = ~paste('P:', p, '; H:', h), hoverinfo = 'text')


  # Add wimp.graphII dots (with normal opacity)
  colorsII <- .construct.colors(wimp = wimp.t1, mode = "red/green")
  phm.matII.df$color <- colorsII[merge$index2,"color"]
  p <- p %>%
    add_markers(data = phm.matII.df, x = ~p, y = ~h,
                marker = list(color = ~color, size = 9,
                              line = list(color = 'black', width = 1)),
                text = ~paste('P:', p, '; H:', h), hoverinfo = 'text')


  # Adding dashed lines between corresponding constructs of wimp.graphI and wimp.graphII
  for (i in 1:nrow(phm.matI.df)) {
    p <- p %>%
      add_segments(x = phm.matI.df$p[i], y = phm.matI.df$h[i],
                   xend = phm.matII.df$p[i], yend = phm.matII.df$h[i],
                   line = list(color = wgI.color, width = 1, dash = 'dash'))
  }

  # Adding construct labels for wimp.graphI
  p <- p %>%
    add_annotations(data = phm.matI.df, x = ~p, y = ~h, text = "",
                    hovertext = ~paste('Constructo:', constructo, '\nP:', p, 'H:', h), hoverinfo = 'text',
                    font = list(size = 12 * text.size, color = wgI.color, opacity = wgI.opacity),
                    showarrow = FALSE, xanchor = 'center', yanchor = 'bottom',
                    yshift = 5)

  # Adding construct labels for wimp.graphII
  p <- p %>%
    add_annotations(data = phm.matII.df, x = ~p, y = ~h, text = ~self.constr,
                    hovertext = ~paste('Constructo:', constructo, '\nP:', p, 'H:', h), hoverinfo = 'text',
                    font = list(size = 12 * text.size, color = 'black'),
                    showarrow = FALSE, xanchor = 'center', yanchor = 'bottom',
                    yshift = 5)

  # Drawing of centroids of both construct systems
  if (show.centroid) {
    centroidI <- phm.matI.df %>%
      summarise(mean_p = mean(p, na.rm = TRUE), mean_h = mean(h, na.rm = TRUE))
    centroidII <- phm.matII.df %>%
      summarise(mean_p = mean(p, na.rm = TRUE), mean_h = mean(h, na.rm = TRUE))

    p <- p %>%
      add_markers(x = centroidI$mean_p, y = centroidI$mean_h,
                  marker = list(color = '#FFD97D', size = 9, symbol = "x", opacity = 0.5,
                                line = list(color = '#FA9D13', width = 1)),
                  text = paste('Centroide Test', 'P:', round(centroidI$mean_p, 5),
                               'H:', format(centroidI$mean_h, nsmall = 5)), hoverinfo = 'text') %>%
      add_markers(x = centroidII$mean_p, y = centroidII$mean_h,
                  marker = list(color = '#FFD97D', size = 12, symbol = "x", opacity = 1,
                                line = list(color = '#FA9D13', width = 1)),
                  text = paste('Centroide Retest', 'P:', round(centroidII$mean_p, 5),
                               'H:', format(centroidII$mean_h, nsmall = 5)), hoverinfo = 'text')
  }

  return(p)
}
