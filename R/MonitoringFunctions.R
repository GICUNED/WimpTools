utils::globalVariables(c("h", "p", ".merge_wimp",
                         ".compatibility_merge_wimp"))

#' Monitoring SSI Adjustment -- monitoring_ssi()
#'
#' @description This function generates an interactive heatmap to visualize
#' and compare the progress of the SSI index across different
#' time points (e.g., pre- and post-intervention).
#'
#' @param wimp_t0 First Subject's WimpGrid (e.g., pre-intervention). It must
#'   be a "wimp" S3 object.
#' @param wimp_t1 Second Subject's WimpGrid (e.g., post-intervention). It must
#'   be a "wimp" S3 object.
#'
#' @return A two heatmap made with plotly.
#'
#' @import plotly
#' @export
#'
#' @examples
#'  monitoring_ssi(example_wimp, example_wimp)

# Monitoring Adjustment Radar Plot ---------------------------------------------

#' Monitoring Self Adjustment -- monitoring_adj()
#'
#' @description This function generates an interactive radar plot to visualize
#' and compare the progress of the "self" and "ideal" across different
#' time points (e.g., pre- and post-intervention).
#'
#' @param wimp_t0 First Subject's WimpGrid (e.g., pre-intervention). It must be
#'   a "wimp" S3 object.
#' @param wimp_t1 Second Subject's WimpGrid (e.g., post-intervention). It must
#'                be a "wimp" S3 object.
#' @param legend If TRUE, displays legend of plot. Default is TRUE.
#'
#' @return Interactive radar plot using plotly
#'
#' @import plotly
#' @export
#'
#' @examples
#'  monitoring_adj(example_wimp, example_wimp)

monitoring_adj <- function(wimp_t0, wimp_t1, legend = TRUE) {

  wimp_t0 <- .align_wimp(wimp_t0, exclude_dilemmatics = FALSE)
  wimp_t1 <- .align_wimp(wimp_t1, exclude_dilemmatics = FALSE)

  merge <- .merge_wimp(wimp_t0, wimp_t1)
  if (.compatibility_merge_wimp(wimp_t0, wimp_t1) == "Incompatibility") {
    stop("WimpGrids have no constructs in common. Monitoring not possible.")
  }

  self_t0 <- wimp_t0$vertices$self[merge$index1]
  self_t0 <- c(self_t0, self_t0[1])

  ideal_t0 <- wimp_t0$vertices$ideal[merge$index1]
  ideal_t0 <- c(ideal_t0, ideal_t0[1])

  right_poles_t0 <- wimp_t0$vertices$right_pole[merge$index1]
  left_poles_t0 <- wimp_t0$vertices$left_pole[merge$index1]
  poles_t0 <- paste(right_poles_t0, " (", left_poles_t0, ")", sep = "")
  poles_t0 <- c(poles_t0, poles_t0[1])

  construct_t0 <- paste(wimp_t0$vertices$left_pole[merge$index1], " - ",
                        wimp_t0$vertices$right_pole[merge$index1], sep = "")
  construct_t0 <- c(construct_t0, construct_t0[1])

  colors_t0 <- .construct_colors(wimp_t0, mode = "red/green")[merge$index1, 1]
  colors_t0 <- c(colors_t0, colors_t0[1])


  self_t1 <- wimp_t1$vertices$self[merge$index2]
  self_t1 <- c(self_t1, self_t1[1])

  ideal_t1 <- wimp_t1$vertices$ideal[merge$index2]
  ideal_t1 <- c(ideal_t1, ideal_t1[1])

  right_poles_t1 <- wimp_t1$vertices$right_pole[merge$index2]
  left_poles_t1 <- wimp_t1$vertices$left_pole[merge$index2]
  poles_t1 <- paste(right_poles_t1, " (", left_poles_t1, ")", sep = "")
  poles_t1 <- c(poles_t1, poles_t1[1])

  construct_t1 <- paste(wimp_t1$vertices$left_pole[merge$index2], " - ",
                        wimp_t1$vertices$right_pole[merge$index2], sep = "")
  construct_t1 <- c(construct_t1, construct_t1[1])

  colors_t1 <- .construct_colors(wimp_t1, mode = "red/green")[merge$index2, 1]
  colors_t1 <- c(colors_t1, colors_t1[1])

  fig <- plot_ly(
    type = "scatterpolar",
    fill = "toself"
  )
  fig <- fig %>%
    add_trace(
      mode = "lines",
      r = 0,
      theta = poles_t0,
      fill = "none",
      line = list(color = "#444444", width = 1.5, shape = "spline",
                  smoothing = 1),
      name = "Pole Threshold",
      hoverinfo = "none",
      showlegend = FALSE
    )
  fig <- fig %>%
    add_trace(
      mode = "lines",
      r = ideal_t0,
      theta = poles_t0,
      fill = "none",
      line = list(color = "darkgreen", width = 2.5, dash = "dot"),
      name = "Baseline Ideal",
      hoverinfo = "none"
    )
  fig <- fig %>%
    add_trace(
      mode = "lines",
      r = ideal_t1,
      theta = poles_t1,
      fill = "none",
      line = list(color = "darkgreen", width = 3),
      name = "Actual Ideal",
      hoverinfo = "none"
    )
  fig <- fig %>%
    add_trace(
      r = self_t0,
      theta = poles_t0,
      name = "Baseline",
      marker = list(color = colors_t0, size = 7,
                    line = list(color = "#FA9D13", width = 1.5)),
      fillcolor = "rgba(255, 217, 125, 0.5)",
      line = list(width = 1, color = "#FA9D13"),
      text = ~paste("<B>", construct_t0, "</B>", "\nSelf:",
                    round(self_t0, 2), "\nIdeal:", round(ideal_t0, 2)),
      hoverinfo = "text",
      hoverlabel = list(bgcolor = colors_t0)
    )

  fig <- fig %>%
    add_trace(
      r = self_t1,
      theta = poles_t1,
      name = "Actual",
      marker = list(color = colors_t1, size = 7,
                    line = list(color = "#AB81E5", width = 1.5)),
      fillcolor = "rgba(213, 192, 242, 0.5)",
      line = list(width = 1, color = "#AB81E5"),
      text = ~paste("<B>", construct_t1, "</B>", "\nSelf:",
                    round(self_t1, 2), "\nIdeal:", round(ideal_t1, 2)),
      hoverinfo = "text",
      hoverlabel = list(bgcolor = colors_t1)
    )
  fig <- fig %>%
    layout(
      showlegend = legend,
      polar = list(
        radialaxis = list(
          visible = TRUE,
          range = if (!is.null(wimp_t0$global$scale)) sort(wimp_t0$global$scale)
          else c(-1, 1)
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
#' @param wimp_t0 First Subject's WimpGrid (e.g., pre-intervention). It must
#'   be a "wimp" S3 object.
#' @param wimp_t1 Second Subject's WimpGrid (e.g., post-intervention). It must
#'   be a "wimp" S3 object.
#'
#' @return A two heatmap made with plotly.
#'
#' @import plotly
#' @export
#'
#' @examples
#'  monitoring_ssi(example_wimp, example_wimp)

monitoring_ssi <- function(wimp_t0, wimp_t1) {

  create_heatmap <- function(wimp, show_y_axis_title = TRUE,
                             show_legend = FALSE, hide_y_ticks = FALSE) {
    x <- wimp$vertices$self
    y <- wimp$vertices$ideal

    alpha_values <- seq(0, 1, by = 0.01)
    beta_values <- seq(0, 1, by = 0.01)

    sim_matrix <- outer(alpha_values, beta_values, Vectorize(
      function(alpha, beta) {
        .sim_index(x, y, alpha = alpha, beta = beta)
      }
    ))

    plot_ly(
      x = alpha_values,
      y = beta_values,
      z = sim_matrix,
      type = "heatmap",
      colorscale = list(
        c(0, "#F52722"),
        c(0.5, "white"),
        c(1, "#A5D610")
      ),
      zmin = 0,
      zmax = 1,
      hovertemplate = paste0(
        "<b>Alpha:</b> %{x}<br>",
        "<b>Beta:</b> %{y}<br>",
        "<b>Adjustment:</b> %{z}<extra></extra>"
      ),
      showscale = show_legend
    ) %>%
      layout(
        xaxis = list(title = "Attention to Self Discrepances"),
        yaxis = list(
          title = if (show_y_axis_title) {
            "Attention to the Desired Change"
          } else {
            NULL
          },
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
            line = list(color = "black", width = 2)
          ),
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
          bgcolor = "rgba(255, 255, 255, 0.8)",
          bordercolor = "black",
          font = list(size = 12)
        )
      )
  }

  heatmap1 <- create_heatmap(wimp_t0, show_y_axis_title = TRUE,
                             show_legend = TRUE) %>%
    layout(title = "Pre-Intervention")
  heatmap2 <- create_heatmap(wimp_t1, show_y_axis_title = FALSE,
                             show_legend = FALSE, hide_y_ticks = TRUE) %>%
    layout(title = "Post-Intervention")

  subplot(heatmap1, heatmap2, nrows = 1, titleX = TRUE, titleY = TRUE,
          margin = 0.005) %>%
    layout(
      title = "",
      showlegend = FALSE
    )
}

# Monitoring PH Index
#' Monitoring PH index -- monitoring_ph()
#'
#' @description This function generates a graphical comparison of constructs
#'   from two different Weighted Implications Grids plotted in a PH space.
#'
#' @param wimp_t0 Data object containing constructs and their respective P and
#'   H coordinates for the first grid.
#' @param wimp_t1 Data object containing constructs and their respective P and
#'   H coordinates for the second grid.
#' @param show_centroid Logical; if TRUE, displays the centroid of construct
#'   P-H coordinates for both grids on the graph.
#' @param text_size Size of the text labels. Default is 1.
#' @param ... additional arguments are passed from \\code{\\link{ph_index}}
#'   function.
#'
#' @return A Plotly object representing the comparative graph of constructs
#'   across both grids.
#'
#' @export
#' @examples
#'  monitoring_ph(example_wimp, example_wimp)

monitoring_ph <- function(wimp_t0, wimp_t1, show_centroid = TRUE,
                          text_size = 1, ...) {

  wimp_t0 <- .align_wimp(wimp_t0, exclude_dilemmatics = FALSE)
  wimp_t1 <- .align_wimp(wimp_t1, exclude_dilemmatics = FALSE)

  merge <- .merge_wimp(wimp_t0, wimp_t1)
  if (.compatibility_merge_wimp(wimp_t0, wimp_t1) == "Incompatibility") {
    stop("WimpGrids have no constructs in common. Monitoring not possible.")
  }

  # Calculate the Presence-Balance index for wimp_t1
  phm_mat_ii <- pb_index(wimp = wimp_t1, ...)[merge$index2, ]
  phm_mat_ii_df <- as.data.frame(phm_mat_ii)
  phm_mat_ii_df$construct <- rownames(phm_mat_ii)
  phm_mat_ii_df$self_constr <- if (!is.null(
    wimp_t1$vertices$self_pole
  )) {
    wimp_t1$vertices$self_pole[merge$index2]
  } else {
    paste(wimp_t1$vertices$left_pole[merge$index2], "/",
          wimp_t1$vertices$right_pole[merge$index2])
  }

  # Calculate the Presence-Balance index for wimp_t0
  phm_mat_i <- pb_index(wimp = wimp_t0, ...)[merge$index1, ]
  phm_mat_i_df <- as.data.frame(phm_mat_i)
  phm_mat_i_df$construct <- rownames(phm_mat_i)
  phm_mat_i_df$self_constr <- if (!is.null(
    wimp_t0$vertices$self_pole
  )) {
    wimp_t0$vertices$self_pole[merge$index1]
  } else {
    paste(wimp_t0$vertices$left_pole[merge$index1], "/",
          wimp_t0$vertices$right_pole[merge$index1])
  }

  # Define the boundaries of the regions
  limit <- max(abs(phm_mat_ii_df$p), abs(phm_mat_ii_df$h)) * 1.1

  # Opacity for wimp_t0 constructs
  wg_i_opacity <- 0.5
  wg_i_color <- "grey"

  # Configuring regions
  shapes <- list(
    list(
      type = "line", x0 = 0, y0 = 0, x1 = limit, y1 = limit,
      xref = "x", yref = "y",
      line = list(color = "#FFD97D", width = 1, dash = "dash")
    ),
    list(
      type = "line", x0 = 0, y0 = 0, x1 = limit, y1 = -limit,
      xref = "x", yref = "y",
      line = list(color = "#FFD97D", width = 1, dash = "dash")
    )
  )

  # Initialise the Plotly chart
  p <- plot_ly()

  # Configuring layout
  p <- p %>%
    layout(
      title = "",
      xaxis = list(title = "Presence"),
      yaxis = list(title = "Hierarchy"),
      plot_bgcolor = "white",
      font = list(family = "Arial"),
      showlegend = FALSE,
      shapes = shapes
    )

  # Add wimp_t0 dots (with lower opacity)
  colors_i <- .construct_colors(wimp = wimp_t0, mode = "red/green")
  phm_mat_i_df$color <- colors_i[merge$index1, "color"]
  p <- p %>%
    add_markers(
      data = phm_mat_i_df, x = ~p, y = ~h,
      marker = list(
        color = ~color, size = 4, opacity = wg_i_opacity,
        line = list(color = "black", width = 1, dash = "dot")
      ),
      text = ~paste("P:", p, "; H:", h), hoverinfo = "text"
    )


  # Add wimp_t1 dots (with normal opacity)
  colors_ii <- .construct_colors(wimp = wimp_t1, mode = "red/green")
  phm_mat_ii_df$color <- colors_ii[merge$index2, "color"]
  p <- p %>%
    add_markers(
      data = phm_mat_ii_df, x = ~p, y = ~h,
      marker = list(
        color = ~color, size = 9,
        line = list(color = "black", width = 1)
      ),
      text = ~paste("P:", p, "; H:", h), hoverinfo = "text"
    )


  # Adding dashed lines between corresponding constructs
  for (i in seq_len(nrow(phm_mat_i_df))) {
    p <- p %>%
      add_segments(
        x = phm_mat_i_df$p[i], y = phm_mat_i_df$h[i],
        xend = phm_mat_ii_df$p[i], yend = phm_mat_ii_df$h[i],
        line = list(color = wg_i_color, width = 1, dash = "dash")
      )
  }

  # Adding construct labels for wimp_t0
  p <- p %>%
    add_annotations(
      data = phm_mat_i_df, x = ~p, y = ~h, text = "",
      hovertext = ~paste("Construct:", construct, "\nP:", p, "H:", h),
      hoverinfo = "text",
      font = list(size = 12 * text_size, color = wg_i_color,
                  opacity = wg_i_opacity),
      showarrow = FALSE, xanchor = "center", yanchor = "bottom",
      yshift = 5
    )

  # Adding construct labels for wimp_t1
  p <- p %>%
    add_annotations(
      data = phm_mat_ii_df, x = ~p, y = ~h, text = ~self_constr,
      hovertext = ~paste("Construct:", construct, "\nP:", p, "H:", h),
      hoverinfo = "text",
      font = list(size = 12 * text_size, color = "black"),
      showarrow = FALSE, xanchor = "center", yanchor = "bottom",
      yshift = 5
    )

  # Drawing of centroids of both construct systems
  if (show_centroid) {
    centroid_i <- phm_mat_i_df %>%
      summarise(mean_p = mean(p, na.rm = TRUE),
                mean_h = mean(h, na.rm = TRUE))
    centroid_ii <- phm_mat_ii_df %>%
      summarise(mean_p = mean(p, na.rm = TRUE),
                mean_h = mean(h, na.rm = TRUE))

    p <- p %>%
      add_markers(
        x = centroid_i$mean_p, y = centroid_i$mean_h,
        marker = list(color = "#FFD97D", size = 9, symbol = "x",
                      opacity = 0.5,
                      line = list(color = "#FA9D13", width = 1)),
        text = paste("Centroid Test", "P:", round(centroid_i$mean_p, 5),
                     "H:", format(centroid_i$mean_h, nsmall = 5)),
        hoverinfo = "text"
      ) %>%
      add_markers(
        x = centroid_ii$mean_p, y = centroid_ii$mean_h,
        marker = list(color = "#FFD97D", size = 12, symbol = "x",
                      opacity = 1,
                      line = list(color = "#FA9D13", width = 1)),
        text = paste("Centroid Retest", "P:", round(centroid_ii$mean_p, 5),
                     "H:", format(centroid_ii$mean_h, nsmall = 5)),
        hoverinfo = "text"
      )
  }

  return(p)
}
