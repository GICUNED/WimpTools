## CHANGE IMPLICATIONS FUNCTIONS

#' Impact and Feedback Index - if_index()
#'
#' @description Computes Impact and Feedback indices for each construct.
#' Impact quantifies influence exerted on other constructs; Feedback
#' quantifies reciprocal influence from the system.
#'
#' @param wimp Subject's WimpGrid. It must be a "wimp" S3 object
#'        imported by the \code{\link{importwimp}} function.
#' @param std Standardization method. One of `"none"`, `"vertex"`,
#'   `"edges"`, or `"adjacent"`. Defaults to `"adjacent"`.
#'
#' @return A data frame with columns:
#'   - "Negative Impact", "Positive Impact", "Global Impact"
#'   - "Negative Feedback", "Positive Feedback", "Global Feedback"
#'   Rows are constructs ("left_pole - right_pole").
#'
#' @details Dilemmatic targets (ideal = 0) are excluded from Impact sums.
#'   Standardization divides results by the chosen coefficient.
#'
#' @author Maite Benitez Santos, Guillermo Calleja Garate,
#'   Alejandro Sanfeliciano
#' @export
#' @examples
#' if_index(example_wimp)
if_index <- function(wimp, std = "adjacent") {
  # Align so all ideal values are positive (required for comparisons)
  wimp <- .align_wimp(wimp, exclude_dilemmatics = FALSE)

  # Base matrices and adjacency (used for standardization when needed)
  wmatrix <- wimp$global$weight_matrix
  amatrix <- wmatrix
  amatrix[amatrix != 0] <- 1

  # Standardization coefficient selection
  if (std == "vertex") std_coef <- nrow(wmatrix) - 1
  if (std == "edges") std_coef <- length(wmatrix[wmatrix != 0])
  if (std == "adjacent") std_coef <- rowSums(amatrix)

  ideal_vector <- wimp$vertices$ideal

  # Impact matrix (zero out dilemmatic targets)
  imatrix <- wmatrix
  dilemmatic_idx <- which(ideal_vector == 0)
  if (length(dilemmatic_idx) > 0) imatrix[, dilemmatic_idx] <- 0

  # Split impact into positive / negative components
  imatrix_pos <- imatrix
  imatrix_pos[imatrix_pos < 0] <- 0
  imatrix_neg <- imatrix
  imatrix_neg[imatrix_pos > 0] <- 0

  negative_impact <- apply(imatrix_neg, 1, sum)
  positive_impact <- apply(imatrix_pos, 1, sum)
  global_impact <- apply(imatrix, 1, sum)

  # Standardize impact if requested
  if (std != "none") {
    std_coef_safe <- std_coef
    std_coef_safe[std_coef_safe == 0] <- 1
    negative_impact <- negative_impact / std_coef_safe
    positive_impact <- positive_impact / std_coef_safe
    global_impact <- global_impact / std_coef_safe
  }

  # Feedback matrix (interaction of bilateral influence)
  fmatrix <- wmatrix * t(wmatrix)
  fmatrix_pos <- fmatrix
  fmatrix_pos[fmatrix_pos < 0] <- 0
  fmatrix_neg <- fmatrix
  fmatrix_neg[fmatrix_pos > 0] <- 0

  negative_feedback <- apply(fmatrix_neg, 1, sum)
  positive_feedback <- apply(fmatrix_pos, 1, sum)
  global_feedback <- apply(fmatrix, 1, sum)

  # Standardize feedback if requested
  if (std != "none") {
    std_coef_safe <- std_coef
    std_coef_safe[std_coef_safe == 0] <- 1
    negative_feedback <- negative_feedback / std_coef_safe
    positive_feedback <- positive_feedback / std_coef_safe
    global_feedback <- global_feedback / std_coef_safe
  }

  # Assemble result data frame
  result <- data.frame(
    abs(negative_impact), positive_impact, global_impact,
    abs(negative_feedback), positive_feedback, global_feedback
  )
  colnames(result) <- c(
    "Negative Impact", "Positive Impact", "Global Impact",
    "Negative Feedback", "Positive Feedback", "Global Feedback"
  )
  constructs <- paste(wimp$vertices$left_pole, "-", wimp$vertices$right_pole)
  rownames(result) <- constructs
  return(result)
}

## IF Plot Function

#' IF Index Scatter Plot - if_plot()
#'
#' @description Scatter plot of Impact (x) vs Feedback (y) per construct.
#' Quadrants highlight typical patterns.
#'
#' @param wimp Subject's WimpGrid. It must be a "wimp" S3 object
#'        imported by the \code{\link{importwimp}} function.
#' @param show Construct filter. One of `"all"`, `"dil"`, `"nodil"`.
#' @param center Axis centering. `"data"` (span data) or `"origin"`.
#' @param text_size Text size multiplier. Defaults to 1.
#' @param ... Additional arguments passed to \code{if_index()}.
#'
#' @return A plotly scatter plot.
#'
#' @details Colors reflect congruency between self and ideal:
#'   green (congruent), red (discrepant), yellow (dilemmatic).
#'
#' @author Maite Benitez Santos, Guillermo Calleja Garate,
#'   Alejandro Sanfeliciano
#' @export
#' @import plotly
#' @examples
#' if_plot(example_wimp)
if_plot <- function(wimp, show = "all", center = "data", text_size = 1, ...) {
  # Align grid and extract core vectors
  wimp <- .align_wimp(wimp, exclude_dilemmatics = FALSE)
  self_vector <- wimp$vertices$self
  ideal_vector <- wimp$vertices$ideal
  left_poles <- wimp$vertices$left_pole
  right_poles <- wimp$vertices$right_pole

  # Derive displayed pole depending on self position
  self_poles <- ifelse(
    self_vector < 0, left_poles,
    ifelse(self_vector > 0, right_poles, paste(left_poles, "-", right_poles))
  )
  constructs <- paste(left_poles, "-", right_poles)
  dil <- which(ideal_vector == 0)

  # Congruency coloring
  congruency <- self_vector / ideal_vector
  construct_color <- sapply(congruency, function(x) {
    if (is.na(x) || is.infinite(x)) {
      "#FFD97D"
    } else if (x < 0) {
      "#d13b43"
    } else if (x > 0) {
      "#5ce75c"
    } else {
      "grey"
    }
  })

  # Build plotting data frame from index results
  df <- if_index(wimp, ...)[c(3, 6)]
  df <- data.frame(df, right_poles, constructs, self_poles, construct_color)
  names(df) <- c("I", "FB", "poles", "construct", "self", "color")
  rownames(df) <- right_poles
  if (show == "nodil") df <- df[-dil, ]
  if (show == "dil") df <- df[dil, ]

  # Axis range selection
  if (center == "data") {
    impact_min <- min(df[1]) - 0.15 * max(abs(df[1]))
    impact_max <- max(df[1]) + 0.15 * max(abs(df[1]))
    feedback_min <- min(df[2]) - 0.15 * max(abs(df[2]))
    feedback_max <- max(df[2]) + 0.15 * max(abs(df[2]))
  }
  if (center == "origin") {
    impact_min <- -(max(abs(df[1])) + 0.15 * max(abs(df[1])))
    impact_max <- max(abs(df[1])) + 0.15 * max(abs(df[1]))
    feedback_min <- -(max(abs(df[2])) + 0.15 * max(abs(df[2])))
    feedback_max <- max(abs(df[2])) + 0.15 * max(abs(df[2]))
  }

  # Quadrant background shapes
  shapes <- list(
    list(
      type = "rect", fillcolor = "palegreen",
      line = list(color = "palegreen"), opacity = 0.3, layer = "below",
      x0 = 0, x1 = 10000, y0 = 0, y1 = 10000
    ),
    list(
      type = "rect", fillcolor = "#eb636b",
      line = list(color = "#eb636b"), opacity = 0.3, layer = "below",
      x0 = 0, x1 = -10000, y0 = 0, y1 = -10000
    ),
    list(
      type = "rect", fillcolor = "#ffe65d",
      line = list(color = "#ffe65d"), opacity = 0.3, layer = "below",
      x0 = 0, x1 = -10000, y0 = 0, y1 = 10000
    ),
    list(
      type = "rect", fillcolor = "#ffe65d",
      line = list(color = "#ffe65d"), opacity = 0.3, layer = "below",
      x0 = 0, x1 = 10000, y0 = 0, y1 = -10000
    )
  )

  # Normalized coordinates for label optimization
  x_range <- range(df$I)
  y_range <- range(df$FB)
  x_span <- diff(x_range)
  y_span <- diff(y_range)
  x_padding <- x_span * 0.1
  y_padding <- y_span * 0.1
  norm_x <- (df$I - (x_range[1] - x_padding)) / (x_span + 2 * x_padding)
  norm_y <- (df$FB - (y_range[1] - y_padding)) / (y_span + 2 * y_padding)
  optimized_positions <- .smart_label_positions(
    x_coords = norm_x, y_coords = norm_y, labels = df$poles,
    distance = 8, text_size = 15 * text_size
  )
  if (nrow(optimized_positions) > 0) {
    optimized_positions$xshift_data <- optimized_positions$xshift *
      (x_span + 2 * x_padding)
    optimized_positions$yshift_data <- optimized_positions$yshift *
      (y_span + 2 * y_padding)
  }

  # Core scatter plot
  fig <- plot_ly(data = df, x = ~I, y = ~FB) %>%
    add_markers(
      data = df, x = ~I, y = ~FB,
      marker = list(
        color = ~color, size = 7,
        line = list(color = "black", width = 1)
      ),
      text = ~ paste(
        "<B>", construct, "</B>", "\nSelf:", self,
        "\nI:", round(I, 2), "\nF:", round(FB, 2)
      ),
      hoverinfo = "text"
    ) %>%
    layout(
      xaxis = list(
        title = "IMPACT", range = c(impact_min, impact_max),
        gridcolor = "white", gridwidth = 0.5, zeroline = TRUE,
        zerolinecolor = "black", zerolinewidth = 2
      ),
      yaxis = list(
        title = "FEEDBACK", range = c(feedback_min, feedback_max),
        gridcolor = "white", gridwidth = 0.5, zeroline = TRUE,
        zerolinecolor = "black", zerolinewidth = 2
      ),
      showlegend = FALSE, shapes = shapes
    )

  # Optimized label annotations
  if (nrow(optimized_positions) > 0) {
    for (i in seq_len(nrow(optimized_positions))) {
      fig <- fig %>% add_annotations(
        x = df$I[i], y = df$FB[i], text = optimized_positions$label[i],
        hoverinfo = "skip", font = list(size = 15 * text_size),
        showarrow = TRUE,
        arrowcolor = "rgba(0,0,0,0.15)",
        arrowwidth = 1,
        arrowsize = 0.5,
        axref = "x",
        ayref = "y",
        ax = optimized_positions$x_data[i],
        ay = optimized_positions$y_data[i],
        xanchor = optimized_positions$xanchor[i],
        yanchor = optimized_positions$yanchor[i]
      )
    }
  }

  fig <- fig %>% .plot_optimization()
  return(fig)
}

## Impact and Feedback Bar Chart Function

#' Impact and Feedback Bar Chart - if_barchart()
#'
#' @description Two-panel horizontal bar chart showing positive and
#' negative components of Impact and Feedback.
#'
#' @param wimp Subject's WimpGrid. It must be a "wimp" S3 object
#'        imported by the \code{\link{importwimp}} function.
#' @param show Construct filter. One of `"all"`, `"dil"`, `"nodil"`.
#' @param ... Additional arguments forwarded to \code{if_index()}.
#'
#' @return A plotly subplot with Impact and Feedback panels.
#'
#' @details Dilemmatic constructs (ideal = 0) are highlighted in yellow
#'   when included.
#'
#' @author Alejandro Sanfeliciano
#' @export
#' @import plotly
#' @examples
#' if_barchart(example_wimp)
if_barchart <- function(wimp, show = "all", lang = "en", ...) {
  # Align and extract base vectors
  wimp <- .align_wimp(wimp, exclude_dilemmatics = FALSE)
  self_vector <- wimp$vertices$self
  ideal_vector <- wimp$vertices$ideal
  left_poles <- wimp$vertices$left_pole
  right_poles <- wimp$vertices$right_pole

  # Display pole based on self orientation
  self_poles <- ifelse(
    self_vector < 0, left_poles,
    ifelse(self_vector > 0, right_poles, paste(left_poles, "-", right_poles))
  )
  constructs <- paste(left_poles, "-", right_poles)
  dil <- which(ideal_vector == 0)

  # Base colors (dilemmatics highlighted in yellow if shown)
  col_p <- rep("#5ce75c", length(constructs))
  col_n <- rep("#d13b43", length(constructs))
  if (show == "all") {
    col_p[dil] <- "#FFD97D"
    col_n[dil] <- "#FFD97D"
  }
  width_line <- rep(1, length(constructs))
  pattern <- rep(0, length(constructs))
  pattern[dil] <- 1

  # Prepare index data (NI/PI/NF/PF + global)
  df <- if_index(wimp, ...)[c(1, 2, 4, 5)]
  df <- data.frame(
    df, df[, 1] + df[, 2], right_poles, self_poles,
    constructs, col_n, col_p, width_line
  )
  df[, 1] <- -df[, 1]
  df[, 3] <- -df[, 3]
  names(df) <- c(
    "NI", "PI", "NF", "PF", "global", "right.poles",
    "self.poles", "construct", "colorn", "colorp",
    "widthline"
  )
  rownames(df) <- right_poles
  if (show == "nodil") df <- df[-dil, ]
  if (show == "dil") df <- df[dil, ]

  # Order df by global impact to simplify JS updates
  df <- df[order(df$global), ]
  df$right.poles <- factor(df$right.poles, levels = df$right.poles)
  
  # Range for symmetric scaling
  range <- max(abs(df[c(1, 2, 3, 4)])) * 1.15

  # Positive / Negative Impact panel
  fig1 <- plot_ly(
    data = df, x = ~PI, y = ~right.poles,
    type = "bar", orientation = "h",
    marker = list(
      color = "#AAF683",
      line = list(
        color = ~colorp,
        width = ~widthline
      ),
      pattern = list(
        shape = ~ ifelse(pattern == 1,
          "/", ""
        ),
        fillmode = "overlay",
        fgcolor = "#FFEE7D",
        size = 20
      )
    ),
    hovertext = ~ paste(
      "<B>", construct, "</B>", "\nSelf:",
      self_poles, "\nPositive Impact:",
      round(PI, 2)
    ), hoverinfo = "text"
  ) %>%
    add_trace(
      x = ~NI, name = "Impact",
      marker = list(
        color = "#EE6055",
        line = list(
          color = ~colorn,
          width = ~widthline
        ),
        pattern = list(
          shape = ~ ifelse(pattern == 1,
            "/", ""
          ),
          fillmode = "overlay",
          fgcolor = "#FFEE7D",
          size = 20
        )
      ),
      hovertext = ~ paste(
        "<B>", construct, "</B>", "\nSelf:",
        self_poles, "\nNegative Impact:",
        round(NI, 2)
      ), hoverinfo = "text"
    ) %>%
    layout(
      barmode = "overlay", bargap = 0.08,
      xaxis = list(
        title = "IMPACT", range = c(-range, range),
        showline = TRUE
      ),
      yaxis = list(title = "", showgrid = TRUE, showline = TRUE),
      showlegend = FALSE
    )

  # Positive / Negative Feedback panel
  fig2 <- plot_ly(
    data = df, x = ~PF, y = ~right.poles,
    type = "bar", orientation = "h", name = "Feedback",
    marker = list(
      color = "#AAF683",
      line = list(
        color = ~colorp,
        width = ~widthline
      ),
      pattern = list(
        shape = ~ ifelse(pattern == 1,
          "/", ""
        ),
        fillmode = "overlay",
        fgcolor = "#FFEE7D",
        size = 20
      )
    ),
    hovertext = ~ paste(
      "<B>", construct, "</B>", "\nSelf:",
      self_poles, "\nPositive Feedback:",
      round(PF, 2)
    ), hoverinfo = "text"
  ) %>%
    add_trace(
      x = ~NF, name = "Feedback",
      marker = list(
        color = "#EE6055",
        line = list(
          color = ~colorn,
          width = ~widthline
        ),
        pattern = list(
          shape = ~ ifelse(pattern == 1,
            "/", ""
          ),
          fillmode = "overlay",
          fgcolor = "#FFEE7D",
          size = 20
        )
      ),
      hovertext = ~ paste(
        "<B>", construct, "</B>", "\nSelf:",
        self_poles, "\nNegative Feedback:",
        round(NF, 2)
      ), hoverinfo = "text"
    ) %>%
    layout(
      barmode = "overlay", bargap = 0.08,
      xaxis = list(
        title = "FEEDBACK", range = c(-range, range),
        showline = TRUE
      ),
      yaxis = list(
        title = "", showgrid = TRUE, showline = TRUE,
        showticklabels = TRUE, side = "right"
      ),
      showlegend = FALSE
    )

  # Combine panels
  fig <- subplot(fig1, fig2, margin = 0.005) %>%
    layout(
      margin = list(r = 100),
      xaxis = list(title = "IMPACT"), yaxis = list(title = ""),
      xaxis2 = list(title = "FEEDBACK"),
      yaxis2 = list(
        title = "", showticklabels = TRUE,
        side = "right", overlaying = "y"
      )
    ) %>%
    .plot_optimization()
    
  t <- wt_i18n(lang)

  # Prepare original arrays for JS in the exact order they are plotted
  if_data <- list(
    dict = t,
    constructs = df$construct,
    right_poles = as.character(df$right.poles),
    self_poles = df$self.poles,
    pi = df$PI,
    ni = df$NI,
    pf = df$PF,
    nf = df$NF,
    pattern = pattern[match(df$right.poles, right_poles)],
    col_rg = .construct_colors(wimp, mode = "red/green")[, "color"][match(df$right.poles, right_poles)],
    col_gs = .construct_colors(wimp, mode = "grey scale")[, "color"][match(df$right.poles, right_poles)],
    col_cb = .construct_colors(wimp, mode = "colorblind")[, "color"][match(df$right.poles, right_poles)],
    col_dk = .construct_colors(wimp, mode = "dark")[, "color"][match(df$right.poles, right_poles)],
    col_pt = .construct_colors(wimp, mode = "pastel")[, "color"][match(df$right.poles, right_poles)],
    col_vd = .construct_colors(wimp, mode = "viridis")[, "color"][match(df$right.poles, right_poles)]
  )

  js_if_panel <- "
    function(el, p_x, data) {
      var x = data;
      var settingsModal = document.createElement('div');
      settingsModal.id = 'if_settings_modal';
      Object.assign(settingsModal.style, {
        position: 'absolute', top: '10px', right: '10px',
        width: '90%', maxWidth: '300px', maxHeight: '80vh', overflowY: 'auto', boxSizing: 'border-box',
        backgroundColor: '#fff', zIndex: '2000', padding: '20px', borderRadius: '8px', 
        boxShadow: '0 4px 20px rgba(0,0,0,0.2)', border: '1px solid #eaeaea', display: 'none', 
        fontFamily: 'Inter, Roboto, sans-serif'
      });
      
      var html = '<div style=\"display:flex; justify-content:space-between; align-items:center; border-bottom:1px solid #eaeaea; padding-bottom:10px; margin-bottom:15px;\">' +
                 '<h3 style=\"margin:0; color:#444; font-size:14px;\">' + (x.dict.vis_options || 'Ajustes') + '</h3>' +
                 '<span id=\"close_if_settings\" style=\"cursor:pointer; font-size:20px; font-weight:bold; color:#888; line-height:1;\">&times;</span>' +
                 '</div>';
                 
      html += '<div style=\"margin-bottom:15px;\">' +
              '<label style=\"display:block; margin-bottom:5px; font-weight:bold; color:#444; font-size:13px;\">' + (x.dict.color_palette || 'Paleta') + '</label>' +
              '<select id=\"if_palette_sel\" style=\"width:100%; padding:4px; border-radius:4px;\">' +
              '<option value=\"rg\" selected>' + (x.dict.pal_redgreen || 'Red-Green') + '</option>' +
              '<option value=\"cb\">' + (x.dict.pal_colorblind || 'Colorblind') + '</option>' +
              '<option value=\"gs\">' + (x.dict.pal_greyscale || 'Greyscale') + '</option>' +
              '<option value=\"dk\">' + (x.dict.pal_dark || 'Dark') + '</option>' +
              '<option value=\"col_pt\">' + (x.dict.pastel || 'Pastel') + '</option>' +
              '<option value=\"col_vd\">' + (x.dict.viridis || 'Viridis') + '</option>' +
              '</select></div>';
              
      html += '<div style=\"margin-bottom:15px;\">' +
              '<label style=\"display:block; margin-bottom:5px; font-weight:bold; color:#444; font-size:13px;\">' + (x.dict.filter_constructs || 'Filtrar Constructos') + '</label>' +
              '<div id=\"if_filter_list\" style=\"max-height:180px; overflow-y:auto; border:1px solid #ddd; padding:5px; border-radius:4px; font-size:12px; background:#f9f9f9;\"></div>' +
              '</div>';
              
      settingsModal.innerHTML = html;
      var container = el.closest('.wt-tab-content') || el.parentElement;
      container.appendChild(settingsModal);
      
      var filterContainer = settingsModal.querySelector('#if_filter_list');
      x.constructs.forEach(function(lbl, idx) {
         var div = document.createElement('div');
         div.style.marginBottom = '4px';
         div.innerHTML = '<label style=\"cursor:pointer; display:flex; align-items:center; color:#555;\"><input type=\"checkbox\" checked value=\"' + idx + '\" class=\"if-construct-cb\" style=\"margin-right:6px; accent-color:#8cc63f;\"> ' + lbl + '</label>';
         filterContainer.appendChild(div);
      });
      
      settingsModal.querySelector('#close_if_settings').onclick = function() { settingsModal.style.display = 'none'; };
      
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
      
      var updateIf = function() {
        var pal = settingsModal.querySelector('#if_palette_sel').value;
        
        var activeIndices = [];
        settingsModal.querySelectorAll('.if-construct-cb').forEach(function(cb) {
           if(cb.checked) activeIndices.push(parseInt(cb.value));
        });
        
        var c_colors = x.col_rg;
        if (pal === 'cb') c_colors = x.col_cb;
        if (pal === 'gs') c_colors = x.col_gs;
        if (pal === 'dk') c_colors = x.col_dk;
        if (pal === 'col_pt') c_colors = x.col_pt;
        if (pal === 'col_vd') c_colors = x.col_vd;
        
        var new_y = [], new_c = [], new_t_pi = [], new_t_ni = [], new_t_pf = [], new_t_nf = [];
        var new_pi = [], new_ni = [], new_pf = [], new_nf = [];
        
        for (var i = 0; i < activeIndices.length; i++) {
           var idx = activeIndices[i];
           new_y.push(x.right_poles[idx]);
           var color = c_colors[idx];
           // Ensure dilemmatics stay yellow if shown, or handle normally
           if (x.pattern[idx] === 1 && pal === 'rg') color = \"#FFD97D\";
           new_c.push(color);
           
           new_pi.push(x.pi[idx]);
           new_ni.push(x.ni[idx]);
           new_pf.push(x.pf[idx]);
           new_nf.push(x.nf[idx]);
           
           new_t_pi.push('<b>' + x.constructs[idx] + '</b><br>Self: ' + x.self_poles[idx] + '<br>Positive Impact: ' + parseFloat(x.pi[idx]).toFixed(2));
           new_t_ni.push('<b>' + x.constructs[idx] + '</b><br>Self: ' + x.self_poles[idx] + '<br>Negative Impact: ' + parseFloat(x.ni[idx]).toFixed(2));
           new_t_pf.push('<b>' + x.constructs[idx] + '</b><br>Self: ' + x.self_poles[idx] + '<br>Positive Feedback: ' + parseFloat(x.pf[idx]).toFixed(2));
           new_t_nf.push('<b>' + x.constructs[idx] + '</b><br>Self: ' + x.self_poles[idx] + '<br>Negative Feedback: ' + parseFloat(x.nf[idx]).toFixed(2));
        }
        
        Plotly.restyle(el, {
           y: [new_y], x: [new_pi], 'marker.line.color': [new_c], hovertext: [new_t_pi]
        }, [0]);
        Plotly.restyle(el, {
           y: [new_y], x: [new_ni], 'marker.line.color': [new_c], hovertext: [new_t_ni]
        }, [1]);
        Plotly.restyle(el, {
           y: [new_y], x: [new_pf], 'marker.line.color': [new_c], hovertext: [new_t_pf]
        }, [2]);
        Plotly.restyle(el, {
           y: [new_y], x: [new_nf], 'marker.line.color': [new_c], hovertext: [new_t_nf]
        }, [3]);
      };
      
      settingsModal.querySelector('#if_palette_sel').onchange = updateIf;
      settingsModal.querySelectorAll('.if-construct-cb').forEach(function(cb) {
         cb.onchange = updateIf;
      });
    }
  "

  fig <- fig %>% htmlwidgets::onRender(js_if_panel, data = if_data)
  return(fig)
}
