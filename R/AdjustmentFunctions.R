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

hypo_plot <- function(wimp, text_size = 1, show_labels = TRUE, lang = "en", ...) {

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
  # Add labels only if requested
  if (show_labels) {
    # Calculate optimal collision-free placements that strictly anchor to dots (no lines)
    layouts <- .calculate_pb_layouts(
      x = df$self, 
      y = df$ideal, 
      labels = df$construct, 
      text_size = text_size
    )
    
    layouts_xshift <- layouts$xshift[[1]]
    layouts_yshift <- layouts$yshift[[1]]
    layouts_xanchor <- layouts$xanchor[[1]]
    layouts_yanchor <- layouts$yanchor[[1]]

    if (nrow(df) > 0) {
      for (i in seq_len(nrow(df))) {
        fig <- fig %>% add_annotations(
          x = df$self[i],
          y = df$ideal[i],
          xshift = layouts_xshift[i],
          yshift = layouts_yshift[i],
          text = df$construct[i],
          showarrow = FALSE,
          xanchor = layouts_xanchor[i],
          yanchor = layouts_yanchor[i],
          font = list(size = 11 * text_size, color = "black"),
          hoverinfo = "skip"
        )
      }
    }
  }

  fig <- fig %>%
    layout(
      margin = list(r = 85),
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
  t <- wt_i18n(lang)

  hypo_data <- list(
    dict = t,
    constructs = df$construct,
    col_rg = .construct_colors(wimp, mode = "red/green")[, "color"],
    col_gs = .construct_colors(wimp, mode = "grey scale")[, "color"],
    col_cb = .construct_colors(wimp, mode = "colorblind")[, "color"],
    col_dk = .construct_colors(wimp, mode = "dark")[, "color"],
    col_pt = .construct_colors(wimp, mode = "pastel")[, "color"],
    col_vd = .construct_colors(wimp, mode = "viridis")[, "color"],
    text_size = text_size,
    orig_x = df$self,
    orig_y = df$ideal,
    layouts_x = list(df$self, df$self, df$self, df$self),
    layouts_y = list(df$ideal, df$ideal, df$ideal, df$ideal)
  )
  
  if (show_labels && exists("layouts_xshift")) {
    hypo_data$layouts_xshift <- layouts$xshift
    hypo_data$layouts_yshift <- layouts$yshift
    hypo_data$layouts_xanchor <- layouts$xanchor
    hypo_data$layouts_yanchor <- layouts$yanchor
  }

  js_hypo_panel <- "
    function(el, p_x, data) {
      var x = data;
      var settingsModal = document.createElement('div');
      settingsModal.id = 'hypo_settings_modal';
      Object.assign(settingsModal.style, {
        position: 'absolute', top: '10px', right: '10px',
        width: '90%', maxWidth: '300px', maxHeight: '80vh', overflowY: 'auto', boxSizing: 'border-box',
        backgroundColor: '#fff', zIndex: '2000', padding: '20px', borderRadius: '8px', 
        boxShadow: '0 4px 20px rgba(0,0,0,0.2)', border: '1px solid #eaeaea', display: 'none', 
        fontFamily: 'Inter, Roboto, sans-serif'
      });
      
      var hypoHTML = '<div style=\"display:flex; justify-content:space-between; align-items:center; border-bottom:1px solid #eaeaea; padding-bottom:10px; margin-bottom:15px;\">' +
                   '<h3 style=\"margin:0; color:#444; font-size:14px;\">' + (x.dict.vis_options || 'Ajustes') + '</h3>' +
                   '<span id=\"close_hypo_settings\" style=\"cursor:pointer; font-size:20px; font-weight:bold; color:#888; line-height:1;\">&times;</span>' +
                   '</div>';
                   
      hypoHTML += '<div style=\"margin-bottom:15px;\">' +
                '<label style=\"display:block; margin-bottom:5px; font-weight:bold; color:#444; font-size:13px;\">' + (x.dict.color_palette || 'Paleta') + '</label>' +
                '<select id=\"hypo_palette_sel\" style=\"width:100%; padding:4px; border-radius:4px;\">' +
                '<option value=\"rg\" selected>' + (x.dict.pal_redgreen || 'Red-Green') + '</option>' +
                '<option value=\"cb\">' + (x.dict.pal_colorblind || 'Colorblind') + '</option>' +
                '<option value=\"gs\">' + (x.dict.pal_greyscale || 'Greyscale') + '</option>' +
                '<option value=\"dk\">' + (x.dict.pal_dark || 'Dark') + '</option>' +
                '<option value=\"col_pt\">' + (x.dict.pastel || 'Pastel') + '</option>' +
                '<option value=\"col_vd\">' + (x.dict.viridis || 'Viridis') + '</option>' +
                '</select></div>';
                
      hypoHTML += '<div style=\"margin-bottom:15px;\">' +
                '<label style=\"display:block; margin-bottom:5px; font-weight:bold; color:#444; font-size:13px;\">' + (x.dict.text_size || 'Tamaño del Texto') + '</label>' +
                '<input type=\"range\" id=\"hypo_text_size\" min=\"0.5\" max=\"2.5\" step=\"0.1\" value=\"' + x.text_size + '\" style=\"width:100%; accent-color:#8cc63f;\">' +
                '</div>';
                
      hypoHTML += '<div style=\"margin-bottom:15px;\">' +
                '<label style=\"display:block; margin-bottom:5px; font-weight:bold; color:#444; font-size:13px;\">' + (x.dict.filter_constructs || 'Filtrar Constructos') + '</label>' +
                '<div id=\"hypo_filter_list\" style=\"max-height:180px; overflow-y:auto; border:1px solid #ddd; padding:5px; border-radius:4px; font-size:12px; background:#f9f9f9;\"></div>' +
                '</div>';
                
      hypoHTML += '<div style=\"margin-bottom:15px;\">' +
                '<label style=\"display:block; margin-bottom:5px; font-weight:bold; color:#444; font-size:13px;\">Etiquetas de Constructos</label>' +
                '<div style=\"margin-bottom:5px; display:flex; gap:10px;\">' +
                '<button id=\"btn_shuffle_labels\" style=\"flex:1; padding:6px; background:#f0f0f0; border:1px solid #ccc; border-radius:4px; cursor:pointer; font-size:12px; transition:0.2s;\">Reordenar</button>' +
                '<button id=\"btn_manual_adj\" style=\"flex:1; padding:6px; background:#f0f0f0; border:1px solid #ccc; border-radius:4px; cursor:pointer; font-size:12px; transition:0.2s;\">Ajuste Manual</button>' +
                '</div></div>';
                
      settingsModal.innerHTML = hypoHTML;
      var container = el.closest('.wt-tab-content') || el.parentElement;
      container.appendChild(settingsModal);
      
      var filterContainer = settingsModal.querySelector('#hypo_filter_list');
      x.constructs.forEach(function(lbl, idx) {
         var div = document.createElement('div');
         div.style.marginBottom = '4px';
         div.innerHTML = '<label style=\"cursor:pointer; display:flex; align-items:center; color:#555;\"><input type=\"checkbox\" checked value=\"' + idx + '\" class=\"hypo-construct-cb\" style=\"margin-right:6px; accent-color:#8cc63f;\"> ' + lbl + '</label>';
         filterContainer.appendChild(div);
      });
      
      settingsModal.querySelector('#close_hypo_settings').onclick = function() { settingsModal.style.display = 'none'; };
      
      var flexbox = el.parentElement.querySelector('div[style*=\"z-index: 1000\"]') || el.parentElement.querySelector('div[style*=\"z-index:1000\"]');
      if (flexbox) {
         var settingsBtn = document.createElement('div');
         settingsBtn.style.cssText = 'background-color:rgba(255,255,255,0.95);width:clamp(26px, 4vmin, 34px);height:clamp(26px, 4vmin, 34px);border-radius:6px;box-shadow:0 2px 10px rgba(0,0,0,0.1);border:1px solid #ddd;display:flex;align-items:center;justify-content:center;cursor:pointer;transition:all 0.2s;';
         settingsBtn.title = x.dict.vis_options || \"Ajustes\";
         settingsBtn.onmouseover = function() { this.style.backgroundColor='#f5f5f5'; };
         settingsBtn.onmouseout = function() { this.style.backgroundColor='rgba(255,255,255,0.95)'; };
         settingsBtn.onclick = function() { settingsModal.style.display = (settingsModal.style.display === 'block' ? 'none' : 'block'); };
         settingsBtn.innerHTML = \"<svg width='60%' height='60%' viewBox='0 0 24 24' fill='none' stroke='#333' stroke-width='2.5' stroke-linecap='round' stroke-linejoin='round'><circle cx='12' cy='12' r='3'></circle><path d='M19.4 15a1.65 1.65 0 0 0 .33 1.82l.06.06a2 2 0 0 1 0 2.83 2 2 0 0 1-2.83 0l-.06-.06a1.65 1.65 0 0 0-1.82-.33 1.65 1.65 0 0 0-1 1.51V21a2 2 0 0 1-2 2 2 2 0 0 1-2-2v-.09A1.65 1.65 0 0 0 9 19.4a1.65 1.65 0 0 0-1.82.33l-.06.06a2 2 0 0 1-2.83 0 2 2 0 0 1 0-2.83l.06-.06a1.65 1.65 0 0 0 .33-1.82 1.65 1.65 0 0 0-1.51-1H3a2 2 0 0 1-2-2 2 2 0 0 1 2-2h.09A1.65 1.65 0 0 0 4.6 9a1.65 1.65 0 0 0-.33-1.82l-.06-.06a2 2 0 0 1 0-2.83 2 2 0 0 1 2.83 0l.06.06a1.65 1.65 0 0 0 1.82.33H9a1.65 1.65 0 0 0 1-1.51V3a2 2 0 0 1 2-2 2 2 0 0 1 2 2v.09a1.65 1.65 0 0 0 1 1.51 1.65 1.65 0 0 0 1.82-.33l.06-.06a2 2 0 0 1 2.83 0 2 2 0 0 1 0 2.83l-.06.06a1.65 1.65 0 0 0-.33 1.82V9a1.65 1.65 0 0 0 1.51 1H21a2 2 0 0 1 2 2 2 2 0 0 1-2 2h-.09a1.65 1.65 0 0 0-1.51 1z'></path></svg>\";
         if (flexbox.children.length > 1) {
            flexbox.insertBefore(settingsBtn, flexbox.children[1]);
         } else {
            flexbox.appendChild(settingsBtn);
         }
      }
      
      var origAnnotations = el.layout.annotations ? JSON.parse(JSON.stringify(el.layout.annotations)) : [];
      var currentLayout = 0;
      
      settingsModal.querySelector('#btn_shuffle_labels').onmouseover = function() { this.style.background='#e4e4e4'; };
      settingsModal.querySelector('#btn_shuffle_labels').onmouseout = function() { this.style.background='#f0f0f0'; };
      settingsModal.querySelector('#btn_shuffle_labels').onclick = function() {
         if (x.layouts_xshift && x.layouts_xshift.length > 0) {
             currentLayout = (currentLayout + 1) % x.layouts_xshift.length;
             updateHypo();
         }
      };
      
      var updateHypo = function() {
        var pal = settingsModal.querySelector('#hypo_palette_sel').value;
        var txtSz = parseFloat(settingsModal.querySelector('#hypo_text_size').value);
        
        var activeIndices = [];
        settingsModal.querySelectorAll('.hypo-construct-cb').forEach(function(cb) {
           if(cb.checked) activeIndices.push(parseInt(cb.value));
        });
        
        var c_colors = x.col_rg;
        if (pal === 'cb') c_colors = x.col_cb;
        if (pal === 'gs') c_colors = x.col_gs;
        if (pal === 'dk') c_colors = x.col_dk;
        if (pal === 'col_pt') c_colors = x.col_pt;
        if (pal === 'col_vd') c_colors = x.col_vd;
        
        var new_x = [], new_y = [], new_c = [], new_t = [];
        for (var i = 0; i < activeIndices.length; i++) {
           var idx = activeIndices[i];
           new_x.push(x.orig_x[idx]);
           new_y.push(x.orig_y[idx]);
           new_c.push(c_colors[idx]);
           new_t.push('<b>' + x.constructs[idx] + '</b><br>Ideal Similarity: ' + parseFloat(x.orig_y[idx]).toFixed(3) + '<br>Self Similarity: ' + parseFloat(x.orig_x[idx]).toFixed(3));
        }
        
        var restyleData = {
           x: [new_x],
           y: [new_y],
           'marker.color': [new_c],
           text: [new_t]
        };
        Plotly.restyle(el, restyleData, [0]);
        
        var newAnnotations = [];
        var seenTexts = {};
        for (var i = 0; i < origAnnotations.length; i++) {
           var ann = JSON.parse(JSON.stringify(origAnnotations[i]));
           if (!ann.text || seenTexts[ann.text]) continue;
           seenTexts[ann.text] = true;
           
           var c_idx = -1;
           for (var k = 0; k < x.constructs.length; k++) {
               if (x.constructs[k] === ann.text) {
                   c_idx = k;
                   break;
               }
           }
           
           if (c_idx !== -1 && x.layouts_xshift && x.layouts_xshift.length > currentLayout) {
               ann.x = x.layouts_x[currentLayout][c_idx];
               ann.y = x.layouts_y[currentLayout][c_idx];
               ann.xshift = x.layouts_xshift[currentLayout][c_idx];
               ann.yshift = x.layouts_yshift[currentLayout][c_idx];
               if (x.layouts_xanchor) ann.xanchor = x.layouts_xanchor[currentLayout][c_idx];
               if (x.layouts_yanchor) ann.yanchor = x.layouts_yanchor[currentLayout][c_idx];
           }
           
           var isActive = false;
           for (var j = 0; j < activeIndices.length; j++) {
               if (x.constructs[activeIndices[j]] === ann.text) {
                   isActive = true;
                   break;
               }
           }
           if (isActive) {
               ann.font.size = 11 * txtSz;
               newAnnotations.push(ann);
           }
        }
        Plotly.relayout(el, { annotations: newAnnotations });
      };
      
      settingsModal.querySelector('#hypo_palette_sel').onchange = updateHypo;
      settingsModal.querySelector('#hypo_text_size').oninput = updateHypo;
      settingsModal.querySelectorAll('.hypo-construct-cb').forEach(function(cb) {
         cb.onchange = updateHypo;
      });
      
      var isManual = false;
      var btnManual = settingsModal.querySelector('#btn_manual_adj');
      btnManual.onmouseover = function() { if(!isManual) this.style.background='#e4e4e4'; };
      btnManual.onmouseout = function() { if(!isManual) this.style.background='#f0f0f0'; };
      btnManual.onclick = function() {
          isManual = !isManual;
          this.style.background = isManual ? '#d0ebd0' : '#e4e4e4';
          this.style.borderColor = isManual ? '#88c588' : '#ccc';
          var currentConfig = (el._fullLayout && el._fullLayout._modeBar) ? el._fullLayout._modeBar.config : (el._context || { displayModeBar: false });
          var newConfig = Object.assign({}, currentConfig, { edits: { annotationPosition: isManual } });
          Plotly.react(el, el.data, el.layout, newConfig);
      };
    }
  "
  
  fig <- fig %>% htmlwidgets::onRender(js_hypo_panel, data = hypo_data)

  return(fig)
}
