## ADJUSTMENT FUNCTIONS ##

utils::globalVariables(c(".estimate_ssi_parameters"))

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
#' @references
#' Sanfeliciano, A., Hurtado-Martínez, C., Botella García del Cid, L., & Saúl, L. A. (2025). Similarity Self/Ideal Index (SSI): A Feature-Based Approach to Modeling Psychological Well-Being. Mathematics.
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
    hypo_val <- .calc_hypo(self_val, ideal_val)
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
  
  # Calculate structural coefficients if using SSI method
  structural_coefs <- NULL
  if (method == "ssi") {
    structural_coefs <- .calc_structural_coefs(wimp)
  }

  result$global <- df_global
  result$construct <- df_construct
  result$structural_coefs <- structural_coefs
  result$method <- c(method, rc_text)
  result$wimp <- wimp
  class(result) <- "self_index"

  return(result)
}

# Self Radar Chart -------------------------------------------------------------

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

  wimp <- .align_wimp(wimp, exclude_dilemmatics = FALSE)

  self <- wimp$vertices$self
  self <- c(self, self[1])

  ideal <- wimp$vertices$ideal
  ideal <- c(ideal, ideal[1])

  r_poles <- wimp$vertices$right_pole
  l_poles <- wimp$vertices$left_pole
  poles <- paste(r_poles, " (", l_poles, ")", sep = "")
  poles <- c(poles, poles[1])

  construct <- .construct_names(wimp)
  construct <- c(construct, construct[1])

  colors <- .construct_colors(wimp, mode = "red/green")[, 1]
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
      name = "Self",
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
      showlegend = TRUE,
      polar = list(
        radialaxis = list(
          visible = TRUE,
          range = c(-1.2, 1),
          tickvals = seq(-1, 1, by = 0.2)
        )
      )
    )

  plot
}

# SSI Heatmap ------------------------------------------------------------------

#'
#' SSI Heatmap -- ssi_heatmap()
#'
#' @description Creates a heatmap showing SSI values across different
#'              alpha (discrepancy salience) and beta (aspiration salience)
#'              parameter combinations. Optionally, estimates alpha and beta
#'              parameters from the data using a probabilistic model based on
#'              construct preference weights.
#'
#' @param wimp Subject's WimpGrid. It must be a "wimp" S3 object
#'        imported by the \code{\link{importwimp}} function.
#' @param estimation Logical. If \code{TRUE}, estimates alpha and beta parameters
#'   using a probabilistic model based on construct preferences (congruence,
#'   discrepancy, and aspiration). If \code{FALSE} (default), uses fixed values
#'   of 0.5 for both parameters.
#' @param palette Character. Color palette for the heatmap. Options: \code{"redgreen"}
#'   (default, red-white-green), \code{"viridis"} (purple to yellow),
#'   \code{"plasma"} (purple to yellow), \code{"inferno"} (black to yellow),
#'   \code{"magma"} (black to white), \code{"cividis"} (blue to yellow),
#'   \code{"turbo"} (blue to red), \code{"picnic"}, \code{"bluered"},
#'   or \code{"rdbu"} (red-white-blue diverging).
#'
#' @return A plotly heatmap.
#'
#' @author Alejandro Sanfeliciano
#'
#' @references
#' Sanfeliciano, A., Hurtado-Martínez, C., Botella García del Cid, L., & Saúl, L. A. (2025). Similarity Self/Ideal Index (SSI): A Feature-Based Approach to Modeling Psychological Well-Being. Mathematics.
#'
#' @import plotly
#' @export
#'
#' @examples
#' ssi_heatmap(example_wimp)
#' ssi_heatmap(example_wimp, estimation = TRUE)
#' ssi_heatmap(example_wimp, palette = "viridis")
#'

ssi_heatmap <- function(wimp, estimation = FALSE, palette = "Redgreen") {
  
  # Extract self and ideal ratings
  x <- wimp$vertices$self
  y <- wimp$vertices$ideal
  
  # Create fixed grid [0,1] for SSI parameter space
  alpha_values <- seq(0, 1, by = 0.01)
  beta_values  <- seq(0, 1, by = 0.01)
  
  # Compute SSI surface across parameter grid
  sim_matrix <- outer(alpha_values, beta_values,
                      Vectorize(function(alpha, beta) {
                        .sim_index(x, y, alpha = alpha, beta = beta)
                      }))
  
  # Estimate probability distribution if requested
  params <- NULL
  pdf_matrix <- NULL
  
  # Always calculate structural coefficients
  structural_coefs <- .calc_structural_coefs(wimp)
  
  if (estimation) {
    params <- .estimate_ssi_parameters(wimp)
    pdf_matrix <- outer(alpha_values, beta_values, params$pdf_function)
  }
  
  # Select color palette

    if (palette == "Redgreen") {
      palette <- list(c(0, "#F52722"), c(0.5, "white"), c(1, "#A5D610"))
    }
  
  # Create base heatmap layer with SSI surface
  plot <- plot_ly() %>%
    add_heatmap(
      x = alpha_values,
      y = beta_values,
      z = t(sim_matrix),
      colorscale = palette,
      zmin = 0,
      zmax = 1,
      hovertemplate = paste("<b>Alpha:</b> %{x}<br><b>Beta:</b> %{y}",
                            "<br><b>SSI:</b> %{z:.3f}<extra></extra>"),
      colorbar = list(title = "<b>SSI</b>")
    )
  
  # Add probability density contours if estimation enabled
  if (estimation && !is.null(pdf_matrix)) {
    
    # Normalize probability density to [0, 1] range
    max_dens <- max(pdf_matrix, na.rm = TRUE)
    pdf_normalized <- pdf_matrix / max_dens
    
    plot <- plot %>%
      add_contour(
        x = alpha_values,
        y = beta_values,
        z = t(pdf_normalized),
        showscale = FALSE,
        contours = list(
          coloring = 'lines',
          start = 0.1,
          end = 0.9,
          size = 0.15,
          showlabels = TRUE,
          labelfont = list(size = 10, color = "black")
        ),
        line = list(color = 'black', width = 2),
        hoverinfo = "skip"
      ) %>%
      add_markers(
        x = params$mu_alpha,
        y = params$mu_beta,
        marker = list(
          color = "#FDE725FF",
          size = 12, 
          line = list(width = 1, color = "black"),
          symbol = "circle"
        ),
        customdata = .sim_index(x, y, alpha = params$mu_alpha, beta = params$mu_beta),
        hovertemplate = paste(
          "<b>Most Likely Profile</b><br>",
          "Alpha: %{x:.3f}<br>",
          "Beta: %{y:.3f}<br>",
          "SSI: %{customdata:.3f}<extra></extra>"
        ),
        showlegend = FALSE
      )
  }
  
  # Configure plot layout and styling
  plot <- plot %>%
    layout(
      title = "",
      xaxis = list(
        title = list(text = "<b>Discrepancy Salience (Alpha)</b>",
                     font = list(size = 18)),
        tickfont = list(size = 14),
        range = c(0, 1),
        constrain = "domain"
      ),
      yaxis = list(
        title = list(text = "<b>Aspiration Salience (Beta)</b>",
                     font = list(size = 18)),
        tickfont = list(size = 14),
        range = c(0, 1),
        scaleanchor = "x",
        scaleratio = 1
      ),
      shapes = list(
        # Plot border
        list(type = "rect", x0 = 0, x1 = 1, y0 = 0, y1 = 1, 
             line = list(color = "black", width = 2)),
        # Diagonal line: Alpha + Beta = 1
        list(type = "line", x0 = 0, y0 = 1, x1 = 1, y1 = 0,
             line = list(color = "black", width = 1, dash = "dot"))
      ),
      annotations = if (!is.null(structural_coefs)) {
        coef_text <- sprintf(
          "<b>ω<sub>α</sub>:</b> %.3f<br><b>ω<sub>β</sub>:</b> %.3f",
          structural_coefs$omega_alpha,
          structural_coefs$omega_beta
        )
        list(
          list(
            x = 0.5, y = 0.5, text = "+", showarrow = FALSE,
            font = list(color = "rgba(0,0,0,1)", size = 25)
          ),
          list(
            x = 0.95,
            y = 0.05,
            text = coef_text,
            showarrow = FALSE,
            xanchor = "right",
            yanchor = "bottom",
            font = list(size = 12, color = "black"),
            bgcolor = "rgba(255, 255, 255, 0.8)",
            bordercolor = "black",
            borderwidth = 1,
            borderpad = 4
          )
        )
      } else {
        list(
          list(
            x = 0.5, y = 0.5, text = "+", showarrow = FALSE,
            font = list(color = "rgba(0,0,0,1)", size = 25)
          )
        )
      }
    ) %>%
    config(displayModeBar = FALSE)

  return(plot)
}

# Hypothetical Scenarios Plot --------------------------------------------------

#' Hypothetical Scenarios Plot  -- hypo_plot()
#'
#' @description Creates a scatter plot showing self-hypothetical similarity
#'              (SHS) vs self-ideal hypothetical similarity (SHI) for each
#'              construct, with congruence color coding.
#'
#' @param wimp Subject's WimpGrid. It must be a "wimp" S3 object
#'        imported by the \code{\link{importwimp}} function.
#' @param text_size Scalar that modifies the text size. Default is 1.
#' @param show_labels Logical. Whether to show construct labels on the plot
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

hypo_plot <- function(wimp, text_size = 1, show_labels = TRUE, ...) {

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
    hypo_val <- .calc_hypo(self_val, ideal_val)
    if (is.na(hypo_val)) {
      pole_names[i] <- paste(wimp$vertices$left_pole[i], "-",
                             wimp$vertices$right_pole[i])
    } else if (hypo_val > 0) {
      pole_names[i] <- wimp$vertices$right_pole[i]
    } else {
      pole_names[i] <- wimp$vertices$left_pole[i]
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
  if (show_labels) {
    # Use smart label positioning to avoid overlaps
    label_positions <- .smart_label_positions(
      x_coords = df$self,
      y_coords = df$ideal,
      labels = df$construct,
      distance = 8,
      text_size = 11 * text_size
    )

    df$xanchor <- label_positions$xanchor
    df$yanchor <- label_positions$yanchor
    df$opt_x <- label_positions$x
    df$opt_y <- label_positions$y

    fig <- fig %>%
      add_annotations(
        data = df,
        x = ~self,
        y = ~ideal,
        text = ~construct,
        hoverinfo = "skip",
        font = list(size = 11 * text_size, color = "black"),
        showarrow = TRUE,
        arrowcolor = "rgba(0,0,0,0.15)",
        arrowwidth = 1,
        arrowsize = 0.5,
        axref = "x",
        ayref = "y",
        ax = ~opt_x,
        ay = ~opt_y,
        xanchor = ~xanchor,
        yanchor = ~yanchor
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
