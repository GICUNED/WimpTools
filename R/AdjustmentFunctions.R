## ADJUSTMENT FUNCTIONS ##

# Construct proportions index -------------------------------------------------

#' Construct Congruence Analysis -- construct_index()
#'
#' @description Calculates frequency and proportion of construct types
#'              (congruent, discrepant, dilemmatic and undefined) based on
#'              self-ideal relationship analysis.
#'
#' @param wimp Subject's WimpGrid. It must be a "wimp" S3 object
#'        imported by the \code{\link{importwimp}} function.
#'
#' @return A matrix with frequency and proportion values for each construct
#'         type: Congruents, Discrepants, Dilemmatics, and Undefined.
#'
#' @author Alejandro Sanfeliciano
#'
#' @export
#'
#' @examples
#' construct_index(example_wimp)
#'

construct_index <- function(wimp) {

  self_vec <- wimp$vertices$self
  ideal_vec <- wimp$vertices$ideal
  congruents <- which(sign(self_vec) == sign(ideal_vec) &
                        self_vec != 0 & ideal_vec != 0)
  discrepants <- which(sign(self_vec) != sign(ideal_vec) &
                         self_vec != 0 & ideal_vec != 0)
  dilemmatics <- which(ideal_vec == 0)
  undefined <- which(self_vec == 0)
  n_congruents <- length(congruents)
  n_discrepants <- length(discrepants)
  n_dilemmatics <- length(dilemmatics)
  n_undefined <- length(undefined)

  n <- nrow(wimp$vertices)

  congruents <- c(n_congruents, n_congruents / n)
  discrepants <- c(n_discrepants, n_discrepants / n)
  dilemmatics <- c(n_dilemmatics, n_dilemmatics / n)
  undefined <- c(n_undefined, n_undefined / n)

  result <- rbind(congruents, discrepants, dilemmatics, undefined)

  rownames(result) <- c("Congruents", "Discrepants", "Dilemmatics", "Undefined")
  colnames(result) <- c("Frequency", "Proportion")

  result
}


# Self Analysis ----------------------------------------------------------------

#' Self-Ideal Analysis -- self_index()
#'
#' @description Calculates global and construct-specific adjustment indices
#'              between self and hypothetical scenarios. Uses SSI Index or
#'              standard correlation methods.
#'
#' @param wimp Subject's WimpGrid. It must be a "wimp" S3 object
#'        imported by the \code{\link{importwimp}} function.
#' @param method Correlation method: "ssi" (SSI Index), "pearson", "kendall",
#'        or "spearman". Default is "ssi".
#' @param rc Use Cohen's rc (reflection invariant). Default is TRUE.
#' @param alpha Discrepancy salience for SSI (0-1). Default is 0.5.
#' @param beta Aspiration salience for SSI (0-1). Default is 0.5.
#'
#' @return List with \code{global} indices (Self/Ideal, Self/Hypo, Ideal/Hypo)
#'         and \code{construct} data frame with congruence classifications
#'         and similarity measures.
#'
#' @author Alejandro Sanfeliciano
#'
#' @export
#'
#' @examples
#' self_index(example_wimp)
#'
#' # Using different methods
#' self_index(example_wimp, method = "pearson")
#' self_index(example_wimp, method = "kendall", rc = FALSE)
#'

self_index <- function(wimp, method = "ssi", rc = TRUE, alpha = .5, beta = .5) {

  result <- list()

  congruence <- character(nrow(wimp$vertices))
  for (i in seq_len(nrow(wimp$vertices))) {
    self_val <- wimp$vertices$self[i]
    ideal_val <- wimp$vertices$ideal[i]
    hypo_val <- .calc.hypo(self_val, ideal_val)
    if (is.na(ideal_val) || ideal_val == 0) {
      congruence[i] <- "Dilemmatic"
    } else if (is.na(hypo_val)) {
      congruence[i] <- "Undefined"
    } else if (sign(ideal_val) == sign(hypo_val)) {
      congruence[i] <- "Congruent"
    } else {
      congruence[i] <- "Discrepant"
    }
  }

  hypo_matrix <- wimp$global$hypo_matrix
  self_vector <- wimp$vertices$self
  ideal_vector <- wimp$vertices$ideal
  hypo_matrix_full <- cbind(self_vector, hypo_matrix, ideal_vector)
  colnames(hypo_matrix_full)[c(1, ncol(hypo_matrix_full))] <- c("SELF", "IDEAL")
  ncol_matrix <- ncol(hypo_matrix_full)
  hypo_names <- colnames(hypo_matrix_full)[-c(1, ncol_matrix)]

  rc_text <- "no rc"
  if (rc) {
    hypo_matrix_full <- rbind(hypo_matrix_full, -hypo_matrix_full)
    rc_text <- "rc"
  }

  if (method == "ssi") {
    self_vector <- hypo_matrix_full[, 1]
    ideal_vector <- hypo_matrix_full[, ncol_matrix]

    self_cor <- apply(hypo_matrix_full[, -c(1, ncol(hypo_matrix_full))], 2,
                      function(col) {
                        .sim_index(self_vector, col,
                                   alpha = alpha, beta = beta)
                      })
    ideal_cor <- apply(hypo_matrix_full[, -c(1, ncol(hypo_matrix_full))], 2,
                       function(col) {
                         .sim_index(ideal_vector, col,
                                    alpha = alpha, beta = beta)
                       })

    self_ideal_cor <- .sim_index(self_vector, ideal_vector,
                                 alpha = alpha, beta = beta)
    self_hypo_cor <- mean(self_cor)
    ideal_hypo_cor <- mean(ideal_cor)

  } else {
    self_vector <- hypo_matrix_full[, 1]
    ideal_vector <- hypo_matrix_full[, ncol_matrix]

    self_cor <- cor(self_vector, hypo_matrix_full[, -c(1, ncol_matrix)],
                    method = method)
    ideal_cor <- cor(ideal_vector, hypo_matrix_full[, -c(1, ncol_matrix)],
                     method = method)

    self_ideal_cor <- cor(self_vector, ideal_vector, method = method)
    self_hypo_cor <- mean(self_cor)
    ideal_hypo_cor <- mean(ideal_cor)
  }

  df_global <- data.frame(
    self_ideal_cor, self_hypo_cor, ideal_hypo_cor
  )
  names(df_global) <- c("Self/Ideal", "Self/Hypo", "Ideal/Hypo")

  df_construct <- data.frame(
    Hypothetical_Scenario = hypo_names,
    Congruence_Scenario = congruence,
    Self_Similarity = round(as.numeric(self_cor), 4),
    Ideal_Similarity = round(as.numeric(ideal_cor), 4),
    stringsAsFactors = FALSE
  )

  names(df_construct) <- c("Hypothetical Scenario", "Congruence Scenario",
                           "SHS", "SHI")

  result$global <- df_global
  result$construct <- df_construct
  result$method <- c(method, rc_text)
  result$wimp <- wimp
  class(result) <- "self_index"

  return(result)
}

# Self Radar Chart ---------------------------------------------------

#' Self-Ideal Radar Chart -- self_plot()
#'
#' @description Creates a radar chart displaying self and ideal ratings
#'              for each construct, showing congruence patterns visually.
#'
#' @param wimp Subject's WimpGrid. It must be a "wimp" S3 object
#'        imported by the \code{\link{importwimp}} function.
#'
#' @return Interactive Plotly radar chart with self (blue) and ideal (green)
#'         traces, including SSI Index in the legend.
#'
#' @author Alejandro Sanfeliciano
#'
#' @import plotly
#' @export
#'
#' @examples
#' self_plot(example_wimp)
#'

self_plot <- function(wimp) {

  wimp <- .align.wimp(wimp, exclude.dilemmatics = FALSE)

  self <- wimp$vertices$self
  self <- c(self, self[1])

  ideal <- wimp$vertices$ideal
  ideal <- c(ideal, ideal[1])

  r_poles <- wimp$vertices$rpole
  l_poles <- wimp$vertices$lpole
  poles <- paste(r_poles, " (", l_poles, ")", sep = "")
  poles <- c(poles, poles[1])

  construct <- paste(wimp$vertices$lpole, "-", wimp$vertices$rpole, sep = " ")
  construct <- c(construct, construct[1])

  colors <- .construct.colors(wimp, mode = "red/green")[, 1]
  colors <- c(colors, colors[1])

  plot <- plot_ly(
    type = "scatterpolar",
    mode = "lines+markers",
    fill = "toself"
  )
  plot <- plot %>%
    add_trace(
      mode = "lines",
      r = 0,
      theta = poles,
      fill = "none",
      line = list(color = "#444444", width = 1.5,
                  shape = "spline", smoothing = 1),
      name = "Pole Threshold",
      hoverinfo = "none"
    )
  plot <- plot %>%
    add_trace(
      mode = "lines",
      r = ideal,
      theta = poles,
      fill = "none",
      line = list(color = "darkgreen", width = 3, shape = "line"),
      name = "Ideal",
      hoverinfo = "none"
    )
  plot <- plot %>%
    add_trace(
      mode = "lines+markers",
      r = self,
      theta = poles,
      name = paste("SSI Index:", round(self_index(wimp)$global[1], 2)),
      marker = list(color = colors, size = 7,
                    line = list(color = "#6F6BFF", width = 1.5)),
      fillcolor = "rgba(204, 203, 248, 0.5)",
      line = list(width = 1, color = "#6F6BFF"),
      text = ~paste("<B>", construct, "</B>", "\nSelf:", round(self, 2),
                    "\nIdeal:", round(ideal, 2)),
      hoverinfo = "text",
      hoverlabel = list(bgcolor = colors)
    )
  plot <- plot %>%
    layout(
      showlegend = FALSE,
      polar = list(
        radialaxis = list(
          visible = TRUE,
          range = c(-1, 1)
        )
      )
    )

  plot
}

# SSI Heatmap -----------------------------------------------------------
#'
#' SSI Heatmap -- ssi_heatmap()
#'
#' @description Creates a heatmap showing SSI values across different
#'              alpha (discrepancy salience) and beta (aspiration salience)
#'              parameter combinations.
#'
#' @param wimp Subject's WimpGrid. It must be a "wimp" S3 object
#'        imported by the \code{\link{importwimp}} function.
#'
#' @return A plotly heatmap.
#'
#' @author Alejandro Sanfeliciano
#'
#' @import plotly
#' @export
#'
#' @examples
#' ssi_heatmap(example_wimp)
#'

ssi_heatmap <- function(wimp) {

  x <- wimp$vertices$self
  y <- wimp$vertices$ideal

  alpha_values <- seq(0, 1, by = 0.01)
  beta_values <- seq(0, 1, by = 0.01)

  sim_matrix <- outer(alpha_values, beta_values,
                      Vectorize(function(alpha, beta) {
                        .sim_index(x, y, alpha = alpha, beta = beta)
                      }))

  plot <- plot_ly(
    x = alpha_values,
    y = beta_values,
    z = t(sim_matrix),
    type = "heatmap",
    colorscale = list(c(0, "#F52722"), c(0.5, "white"), c(1, "#A5D610")),
    zmin = 0,
    zmax = 1,
    hovertemplate = paste("<b>Alpha:</b> %{x}<br><b>Beta:</b> %{y}",
                          "<br><b>Similarity:</b> %{z}<extra></extra>"),
    colorbar = list(
      title = "<b>SSI</b>",
      tickfont = list(size = 16),
      ticklen = 10
    )
  ) %>%
    layout(
      title = "",
      xaxis = list(
        title = list(
          text = "<b>Discrepancy Salience (Alpha)</b>",
          font = list(size = 25)
        ),
        tickfont = list(size = 18)
      ),
      yaxis = list(
        title = list(
          text = "<b>Aspiration Salience (Beta)</b>",
          font = list(size = 25)
        ),
        tickfont = list(size = 18)
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

  plot
}

#' Hypothetical Scenarios Plot  -- hypo_plot()
#'
#' @description Creates a scatter plot showing self-hypothetical similarity
#'              (SHS) vs self-ideal hypothetical similarity (SHI) for each
#'              construct, with congruence color coding.
#'
#' @param wimp Subject's WimpGrid. It must be a "wimp" S3 object
#'        imported by the \code{\link{importwimp}} function.
#' @param text.size Scalar that modifies the text size. Default is 1.
#' @param show.labels Logical. Whether to show construct labels on the plot
#'        Default is TRUE. Set to FALSE to reduce visual clutter with many
#'        constructs.
#' @param ... Additional arguments passed to \code{\link{self_index}} function.
#'
#'
#' @author  Alejandro Sanfeliciano
#'
#' @return returns a interactive scatter plot made with Plotly.
#'
#' @export
#'
#' @import plotly
#'
#' @examples
#' hypo_plot(example_wimp)
#'
#' # Without labels for cleaner view
#' hypo_plot(example_wimp, show.labels = FALSE)

hypo_plot <- function(wimp, text.size = 1, show.labels = TRUE, ...) { 

  hypo_matrix <- wimp$global$hypo_matrix
  self_vector <- wimp$vertices$self
  ideal_vector <- wimp$vertices$ideal
  hypo_matrix_full <- cbind(self_vector, hypo_matrix, ideal_vector)
  colnames(hypo_matrix_full)[c(1, ncol(hypo_matrix_full))] <- c("SELF", "IDEAL")

  self_index_data <- self_index(wimp, ...)

  congruence <- self_index_data$construct[[2]]

  construct_color <- ifelse(
    congruence == "Congruent", "#A5D610",
    ifelse(congruence == "Discrepant", "#F52722",
           ifelse(congruence == "Undefined", "yellow",
                  ifelse(congruence == "Dilemmatic", "yellow", "#000000")))
  )

  pole_names <- character(nrow(wimp$vertices))
  for (i in seq_len(nrow(wimp$vertices))) {
    self_val <- wimp$vertices$self[i]
    ideal_val <- wimp$vertices$ideal[i]
    hypo_val <- .calc.hypo(self_val, ideal_val)
    if (is.na(hypo_val)) {
      pole_names[i] <- paste(wimp$vertices$lpole[i], "-",
                             wimp$vertices$rpole[i])
    } else if (hypo_val > 0) {
      pole_names[i] <- wimp$vertices$rpole[i]
    } else {
      pole_names[i] <- wimp$vertices$lpole[i]
    }
  }
  df <- self_index_data$construct[c(4, 3)]
  df <- data.frame(df, construct_color, pole_names)
  names(df) <- c("ideal", "self", "color", "construct")
  rownames(df) <- pole_names

  y_ref <- self_index_data[[1]][[1]]

  fig <- plot_ly(
    data = df,
    x = ~self,
    y = ~ideal
  ) %>%
    add_markers(
      data = df,
      x = ~self,
      y = ~ideal,
      marker = list(color = ~color, size = 8,
                    line = list(color = "black", width = 1)),
      text = ~paste("<b>", construct, "</b>",
                    "<br>Ideal Similarity:", round(ideal, 3),
                    "<br>Self Similarity:", round(self, 3)),
      hoverinfo = "text",
      showlegend = FALSE
    )

  # Add labels only if requested
  if (show.labels) {
    # Use smart label positioning to avoid overlaps
    label_positions <- .smart_label_positions(
      x_coords = df$self,
      y_coords = df$ideal,
      labels = df$construct,
      distance = 8,
      text_size = 11 * text.size
    )

    df$xanchor <- label_positions$xanchor
    df$yanchor <- label_positions$yanchor
    df$xshift <- label_positions$xshift
    df$yshift <- label_positions$yshift

    fig <- fig %>%
      add_annotations(
        data = df,
        x = ~self,
        y = ~ideal,
        text = ~construct,
        hoverinfo = "skip",
        font = list(size = 11 * text.size, color = "black"),
        showarrow = FALSE,
        xanchor = ~xanchor,
        xshift = ~xshift,
        yanchor = ~yanchor,
        yshift = ~yshift
      )
  }

  fig <- fig %>%
    layout(
      xaxis = list(
        title = "SELF SIMILARITY (SHS)",
        range = c(0, 1),
        gridwidth = 0.5,
        zeroline = TRUE,
        zerolinecolor = "black",
        zerolinewidth = 2
      ),
      yaxis = list(
        title = "IDEAL SIMILARITY (SHI)",
        range = c(0, 1),
        gridwidth = 0.5,
        zeroline = TRUE,
        zerolinecolor = "black",
        zerolinewidth = 2
      ),
      showlegend = FALSE,
      shapes = list(
        list(
          type = "rect",
          x0 = 0, x1 = 1,
          y0 = y_ref, y1 = 1,
          xref = "x", yref = "y",
          fillcolor = "rgba(0, 255, 0, 0.1)",
          line = list(width = 0)
        ),
        list(
          type = "rect",
          x0 = 0, x1 = 1,
          y0 = 0, y1 = y_ref,
          xref = "x", yref = "y",
          fillcolor = "rgba(255, 0, 0, 0.1)",
          line = list(width = 0)
        ),
        list(
          type = "line",
          x0 = 0, x1 = 1,
          y0 = y_ref, y1 = y_ref,
          xref = "x", yref = "y",
          line = list(color = "#A5D610", dash = "dash", width = 2)
        )
      )
    )

  return(fig)
}
