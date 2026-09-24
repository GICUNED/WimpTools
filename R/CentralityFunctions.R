## CENTRALITY FUNCTIONS ##

# Degree Index Centrality ------------------------------------------------------

#' Degree Centrality Index -- degree_index()
#'
#' @description Calculates centrality based on construct connections. Centrality
#'              represents the degree of connection each construct maintains
#'              with others (number of links per vertex).
#'
#' @param wimp Subject's WimpGrid. It must be a "wimp" S3 object
#'        imported by the \code{\link{importwimp}} function.
#' @param method Centrality calculation method: "simple", "norm", "weight",
#'        "wnorm", or "ego". Default is "weight".
#'
#' @author Alejandro Sanfeliciano
#'
#' @return Returns a list with the centrality data by construct and separated by
#'         input degree, output degree and total degree (in and out).
#'
#' @export
#'
#' @examples
#'
#' degree_index(example_wimp)
#'

degree_index <- function(wimp, method = "weight") {

  poles <- .construct_names(wimp)
  wmat <- wimp$global$weight_matrix

  if (is.null(wmat)) {
    stop("No weights matrix found in wimp object")
  }

  n <- nrow(wmat)

  if (method %in% c("simple", "norm", "ego")) {
    wmat_binary <- wmat / wmat
    wmat_binary[is.nan(wmat_binary)] <- 0
    k_out <- rowSums(wmat_binary)
    k_in <- colSums(wmat_binary)
  }

  if (method %in% c("weight", "wnorm")) {
    k_out <- rowSums(abs(wmat))
    k_in <- colSums(abs(wmat))
  }

  if (method %in% c("norm", "wnorm")) {
    k_out <- k_out / (n - 1)
    k_in <- k_in / (n - 1)
  }

  if (method == "ego") {
    k_out <- k_out / (n * (n - 1))
    k_in <- k_in / (n * (n - 1))
  }

  result <- cbind(k_out, k_in, k_out + k_in)
  rownames(result) <- poles
  colnames(result) <- c("Out", "In", "All")

  return(result)
}

# Distance Matrix --------------------------------------------------------------

#' Shortest Distance Matrix -- dismatrix()
#'
#' @description Calculates shortest distances between construct pairs in the
#'              implication digraph.
#'
#' @param wimp Subject's WimpGrid. It must be a "wimp" S3 object
#'        imported by the \code{\link{importwimp}} function.
#' @param mode Distance calculation mode: "out" (respecting edge direction),
#'        "in" (inverse direction), "all" (ignoring direction). Default is
#'        "out".
#'
#' @author Alejandro Sanfeliciano
#'
#' @return Returns the digraph distance matrix. Matrix that contains the
#'         distances of the shortest paths from one construct to another.
#'
#' @export
#'

dismatrix <- function(wimp, mode = "out") {

  poles <- .construct_names(wimp)
  wmat <- wimp$global$weight_matrix

  if (is.null(wmat)) {
    stop("No weights matrix found in wimp object")
  }

  g <- igraph::graph.adjacency(wmat, mode = "directed", weighted = TRUE)
  result <- igraph::shortest.paths(g, weights = NA, mode = mode)

  rownames(result) <- poles
  colnames(result) <- poles

  return(result)
}

# Closeness Centrality Index ---------------------------------------------------

#' Closeness Centrality Index -- close_index()
#'
#' @description Calculates closeness centrality of constructs within the
#'              implication digraph (inverse of average shortest distance).
#'
#' @param wimp Subject's WimpGrid. It must be a "wimp" S3 object
#'        imported by the \code{\link{importwimp}} function.
#' @param norm If TRUE, values will be normalized. Default is TRUE.
#'
#' @author Alejandro Sanfeliciano
#'
#' @return Returns a vector with the closeness index for each of the
#'         constructs.
#'
#' @export
#'
#' @examples
#'
#' close_index(example_wimp)
#' close_index(example_wimp, norm = FALSE)
#'

close_index <- function(wimp, norm = TRUE) {

  poles <- .construct_names(wimp)
  dist <- dismatrix(wimp)
  n <- nrow(dist)

  if (norm) {
    result <- (n - 1) / rowSums(dist)
  } else {
    result <- 1 / rowSums(dist)
  }

  result <- matrix(result)
  rownames(result) <- poles
  colnames(result) <- "Closeness"

  return(result)
}

# Betweenness Centrality Index -------------------------------------------------

#' Betweenness Centrality Index -- betw_index()
#'
#' @description Calculates betweenness centrality (number of shortest paths
#'              passing through each construct).
#'
#' @param wimp Subject's WimpGrid. It must be a "wimp" S3 object
#'        imported by the \code{\link{importwimp}} function.
#' @param norm If TRUE, values will be normalized. Default is TRUE.
#'
#' @author Alejandro Sanfeliciano
#'
#' @return Returns a vector with the betweeness index for each of the
#'         constructs.
#'
#' @export
#'
#' @examples
#'
#' betw_index(example_wimp)
#' betw_index(example_wimp, norm = FALSE)
#'

betw_index <- function(wimp, norm = TRUE) {

  poles <- .construct_names(wimp)
  wmat <- wimp$global$weight_matrix

  if (is.null(wmat)) {
    stop("No weights matrix found in wimp object")
  }

  g <- igraph::graph.adjacency(wmat, mode = "directed", weighted = TRUE)
  result <- igraph::betweenness(g, normalized = norm, weights = NA)

  result <- matrix(result)
  rownames(result) <- poles
  colnames(result) <- "Betweenness"

  return(result)
}

# PB Centrality Index ----------------------------------------------------------

#' Presence and Balance Indices -- pb_index()
#'
#' @description Computes presence (P) and balance (B) indices for constructs.
#'              P represents frequency of occurrence, B represents influence
#'              balance.
#'
#' @param wimp A wimp object containing implication grid and constructs.
#' @param method Method for calculating degree indices. Default is "wnorm".
#'        Options: "wnorm", "simple", "weight", etc.
#' @param std Standardization method for P and B indices. Options:
#'        \itemize{
#'          \item 'none': No standardization (default)
#'          \item 'vertices': Standardize by max total degree
#'          \item 'edges': Standardize by total number of edges
#'          \item 'max_edges': Standardize by max outgoing edges
#'          \item 'density': Adjust by grid density
#'        }
#'
#' @author Carlos Hurtado and Alejandro Sanfeliciano
#'
#' @return A matrix with two columns, 'p' for presence and 'b' for balance,
#'         containing the indices for each construct. If standardization is
#'         applied, these values are modified according to the selected method.
#'
#' @references
#' Sanfeliciano, A., Saúl, L. A., Hurtado-Martínez, C., & Botella, L. (2025). PB Space: A Mathematical Framework for Modeling Presence and Implication Balance in Psychological Change Through Fuzzy Cognitive Maps. Axioms.
#'
#' @export
#'
#' @examples
#'
#' pb_index(example_wimp)
#' pb_index(example_wimp, std = "vertices")
#' pb_index(example_wimp, method = "wnorm", std = "none")
#'

pb_index <- function(wimp, method = "wnorm", std = "none") {

  c_io <- degree_index(wimp, method = method)
  c_io <- c_io[, c(2, 1, 3)]  # Rearrange In-Out columns
  in_out <- c_io[, 1:2]

  # Linear transformation matrix
  coef <- 1 / sqrt(2)
  coef_matrix <- matrix(c(coef, -coef, coef, coef), nrow = 2)

  # Calculate P-B matrix
  pb_mat <- in_out %*% t(coef_matrix)
  colnames(pb_mat) <- c("p", "b")

  # Standardization
  if (std == "vertices") {
    vertices <- nrow(wimp$vertices)
    coef_max_p <- 2 * coef * (vertices - 1)
    coef_max_b <- coef * (vertices - 1)
    pb_mat[, "p"] <- pb_mat[, "p"] / coef_max_p
    pb_mat[, "b"] <- pb_mat[, "b"] / coef_max_b

  } else if (std == "edges") {
    c_direct_io <- degree_index(wimp, method = "simple")
    c_direct_io <- c_direct_io[, c(2, 1, 3)]
    edges <- sum(c_direct_io[, 2])
    coef_max <- edges * coef
    pb_mat[, "p"] <- pb_mat[, "p"] / coef_max
    pb_mat[, "b"] <- pb_mat[, "b"] / coef_max

  } else if (std == "max_edges") {
    edges <- max(c_io[, 2])
    pb_mat[, "p"] <- pb_mat[, "p"] / edges
    pb_mat[, "b"] <- pb_mat[, "b"] / edges

  } else if (std == "density") {
    vertices <- nrow(wimp$vertices)
    max_edges <- vertices * (vertices - 1)
    c_direct_io <- degree_index(wimp, method = "simple")
    c_direct_io <- c_direct_io[, c(2, 1, 3)]
    total_edges <- sum(c_direct_io[, 2])
    dens <- total_edges / max_edges
    pb_mat[, "p"] <- pb_mat[, "p"] * dens
    pb_mat[, "b"] <- pb_mat[, "b"] * dens
  }

  return(pb_mat)
}

# Eigen Indices ----------------------------------------------------------------

#' Eigenvalue Centrality Index -- eigen_index()
#'
#' @description Calculates centrality scores based on eigenvalue decomposition
#'              of the adjacency matrix.
#'
#' @param wimp An object of class 'wimp' (weighted implications grid).
#' @param matrix Matrix type for analysis: 'direct', 'weights', or
#'        'implications'. Default is 'weights'. Note: only 'weights' is
#'        available in new format.
#' @param num_vectors Number of eigenvectors to use for centrality computation.
#'
#' @author Carlos Hurtado
#'
#' @return A dataframe containing the constructs' names and their respective
#'         centrality scores.
#'
#' @export
#'
#' @examples
#'
#' eigen_index(example_wimp)
#'

eigen_index <- function(wimp, matrix = "weights", num_vectors = 2) {

  if (!matrix %in% c("direct", "weights", "implications")) {
    stop("matrix debe ser 'direct', 'weights' o 'implications'.")
  }

  if (matrix != "weights") {
    warning("Only 'weights' matrix available in new wimp format; using 
            weights.")
  }

  adj_matrix <- wimp$global$weight_matrix

  if (is.null(adj_matrix)) {
    stop("No weights matrix found in wimp object")
  }

  results <- eigen(adj_matrix)

  if (num_vectors > length(results$values)) {
    stop("num.vectors exceeds the number of available eigenvectors.")
  }

  # Calculate centrality using specified number of eigenvectors
  centrality <- Reduce(`+`, lapply(1:num_vectors, function(i) {
    Re(results$vectors[, i])^2 * Re(results$values[i])
  }))

  df_centrality <- data.frame(
    Constructs = .construct_names(wimp),
    Eigenvalues = abs(centrality)
  )

  return(df_centrality)
}

# PB Plot ----------------------------------------------------------------------


.calculate_pb_layouts <- function(x, y, labels, text_size) {
  n <- length(x)
  layouts_xshift <- list()
  layouts_yshift <- list()
  layouts_xanchor <- list()
  layouts_yanchor <- list()
  
  char_width <- 6.5 * text_size
  char_height <- 14 * text_size
  gap <- 1
  
  candidates <- list(
    list(dx = gap+3, dy = 0, xa = "left", ya = "middle"),      # Right (Orthogonal)
    list(dx = -gap-3, dy = 0, xa = "right", ya = "middle"),    # Left (Orthogonal)
    list(dx = 0, dy = gap+2, xa = "center", ya = "bottom"),    # Top (Orthogonal)
    list(dx = 0, dy = -gap-2, xa = "center", ya = "top"),      # Bottom (Orthogonal)
    list(dx = gap, dy = gap, xa = "left", ya = "bottom"),      # Top-Right (Diagonal)
    list(dx = -gap, dy = gap, xa = "right", ya = "bottom"),    # Top-Left (Diagonal)
    list(dx = gap, dy = -gap, xa = "left", ya = "top"),        # Bottom-Right (Diagonal)
    list(dx = -gap, dy = -gap, xa = "right", ya = "top")       # Bottom-Left (Diagonal)
  )
  
  x_range <- max(x) - min(x)
  y_range <- max(y) - min(y)
  if(x_range == 0) x_range <- 1
  if(y_range == 0) y_range <- 1
  
  virt_x <- (x - min(x)) / x_range * 800
  virt_y <- (y - min(y)) / y_range * 600
  
  for (layout_idx in 1:4) {
    xshift <- numeric(n)
    yshift <- numeric(n)
    xanchor <- character(n)
    yanchor <- character(n)
    
    placed_boxes <- list() 
    
    set.seed(layout_idx * 42)
    process_order <- sample(1:n)
    cand_order <- c(layout_idx:8, 1:(layout_idx-1))
    if (layout_idx == 1) cand_order <- 1:8
    
    for (i in process_order) {
      best_cand <- NULL
      best_overlap <- Inf
      
      tw <- nchar(labels[i]) * char_width
      th <- char_height
      
      for (c_idx in cand_order) {
        cand <- candidates[[c_idx]]
        
        left <- virt_x[i] + cand$dx
        if (cand$xa == "right") left <- virt_x[i] + cand$dx - tw
        if (cand$xa == "center") left <- virt_x[i] + cand$dx - tw/2
        
        bottom <- virt_y[i] + cand$dy
        if (cand$ya == "top") bottom <- virt_y[i] + cand$dy - th
        if (cand$ya == "middle") bottom <- virt_y[i] + cand$dy - th/2
        
        right <- left + tw
        top <- bottom + th
        
        overlap <- 0
        for (j in 1:n) {
          if (i == j) next
          if (virt_x[j] >= left && virt_x[j] <= right && virt_y[j] >= bottom && virt_y[j] <= top) {
            overlap <- overlap + 2000
          }
        }
        
        for (box in placed_boxes) {
          if (!(left >= box$right || right <= box$left || bottom >= box$top || top <= box$bottom)) {
             ix_left <- max(left, box$left)
             ix_right <- min(right, box$right)
             iy_bottom <- max(bottom, box$bottom)
             iy_top <- min(top, box$top)
             area <- max(0, ix_right - ix_left) * max(0, iy_top - iy_bottom)
             overlap <- overlap + area
          }
        }
        
        if (overlap < best_overlap) {
          best_overlap <- overlap
          best_cand <- cand
        }
        
        if (overlap == 0) break
      }
      
      xshift[i] <- best_cand$dx
      yshift[i] <- best_cand$dy
      xanchor[i] <- best_cand$xa
      yanchor[i] <- best_cand$ya
      
      left <- virt_x[i] + best_cand$dx
      if (best_cand$xa == "right") left <- left - tw
      if (best_cand$xa == "center") left <- left - tw/2
      bottom <- virt_y[i] + best_cand$dy
      if (best_cand$ya == "top") bottom <- bottom - th
      if (best_cand$ya == "middle") bottom <- bottom - th/2
      
      placed_boxes[[length(placed_boxes) + 1]] <- list(left=left, right=left+tw, bottom=bottom, top=bottom+th)
    }
    
    layouts_xshift[[layout_idx]] <- xshift
    layouts_yshift[[layout_idx]] <- yshift
    layouts_xanchor[[layout_idx]] <- xanchor
    layouts_yanchor[[layout_idx]] <- yanchor
  }
  
  return(list(
    xshift = layouts_xshift,
    yshift = layouts_yshift,
    xanchor = layouts_xanchor,
    yanchor = layouts_yanchor
  ))
}

#' PB Space Scatter Plot -- pb_plot()
#'
#' @description Creates a scatter plot of constructs in Presence-Balance
#'              space. P represents construct frequency, B represents
#'              construct influence balance.
#'
#' @param wimp Subject's WimpGrid. It must be a "wimp" S3 object
#'        imported by the \\code{\\link{importwimp}} function.
#' @param text_size Text label size. Default is 1.
#' @param ... Additional arguments passed to \\code{\\link{pb_index}} function.
#'
#' @author Carlos Hurtado and Alejandro Sanfeliciano
#'
#' @return A Plotly object representing the generated scatter plot.
#'
#' @references
#' Sanfeliciano, A., Saúl, L. A., Hurtado-Martínez, C., & Botella, L. (2025). PB Space: A Mathematical Framework for Modeling Presence and Implication Balance in Psychological Change Through Fuzzy Cognitive Maps. Axioms.
#'
#' @import plotly
#' @export
#'
#' @examples
#'
#' pb_plot(example_wimp)
#'

pb_plot <- function(wimp, text_size = 1, lang = "en", ...) {

  t <- wt_i18n(lang)

  pb_mat <- pb_index(wimp, ...)
  pb_mat_df <- as.data.frame(pb_mat)
  pb_mat_df$construct <- rownames(pb_mat)

  # Get construct labels for display
  if ("self_pole" %in% names(wimp$vertices)) {
    pb_mat_df$self_constr <- wimp$vertices$self_pole
  } else {
    pb_mat_df$self_constr <- .construct_names(wimp)
  }

  # Set plot limits with margin
  limit <- max(abs(pb_mat_df$p), abs(pb_mat_df$b)) * 1.1
  # Presence (p) is a degree-based measure and is never negative, so the
  # x-axis only needs to reach from 0 out to the data (plus a small right
  # margin) - a symmetric autorange left a sliver of unused negative space
  # that became proportionally more visible as the widget got narrower.
  x_right <- max(pb_mat_df$p, na.rm = TRUE) * 1.3
  if (!is.finite(x_right) || x_right <= 0) x_right <- 0.1
  # Round values for display
  pb_mat_df$p <- round(pb_mat_df$p, 3)
  pb_mat_df$b <- round(pb_mat_df$b, 3)

  # Get construct colors
  pb_mat_df$color <- .construct_colors(wimp = wimp,
                                       mode = "red/green")[, "color"]

  # Define plot shapes for non-viable area
  shapes <- list(
    list(type = "path",
         path = paste("M 0,0 L", limit, ",", limit, " L0,", limit, " Z"),
         fillcolor = "#CCCBF8", opacity = 0.2,
         line = list(color = "#CCCBF8")),
    list(type = "path",
         path = paste("M 0,0 L", limit, ",", -limit, " L0,", -limit, " Z"),
         fillcolor = "#CCCBF8", opacity = 0.2,
         line = list(color = "#CCCBF8")),
    list(type = "line", x0 = 0, y0 = 0, x1 = limit, y1 = limit,
         xref = "x", yref = "y",
         line = list(color = "#6F6BFF", width = 1, dash = "dash")),
    list(type = "line", x0 = 0, y0 = 0, x1 = limit, y1 = -limit,
         xref = "x", yref = "y",
         line = list(color = "#6F6BFF", width = 1, dash = "dash"))
  )

  # Calculate optimal collision-free placements that strictly anchor to dots
  layouts <- .calculate_pb_layouts(
    x = pb_mat_df$p, 
    y = pb_mat_df$b, 
    labels = pb_mat_df$self_constr, 
    text_size = text_size
  )
  
  layouts_xshift <- layouts$xshift
  layouts_yshift <- layouts$yshift
  layouts_xanchor <- layouts$xanchor
  layouts_yanchor <- layouts$yanchor
  

  # Create plotly graph
  p <- plot_ly() %>%
    layout(
      title = "",
      xaxis = list(title = list(text = t$pb_axis_x, font = list(size = 20),
                                standoff = 4), range = c(0, x_right),
                  tickformat = ".1f", tickfont = list(size = 13)),
      yaxis = list(title = list(text = t$pb_axis_y,
                                font = list(size = 20), standoff = 4),
                  tickformat = ".1f", tickfont = list(size = 13)),
      margin = list(l = 55, r = 20, t = 10, b = 40),
      plot_bgcolor = "white",
      font = list(family = "Arial"),
      showlegend = FALSE,
      shapes = shapes
    ) %>%
    add_markers(
      data = pb_mat_df, x = ~p, y = ~b,
      mode = "markers",
      marker = list(color = pb_mat_df$color, size = 7,
                    line = list(color = "black", width = 1)),
      text = ~paste("Construct:", self_constr, "<br>P:", p, "<br>B:", b),
      hoverinfo = "text"
    )

  # Add annotations using the first layout
  layouts_x <- list(pb_mat_df$p, pb_mat_df$p, pb_mat_df$p, pb_mat_df$p)
  layouts_y <- list(pb_mat_df$b, pb_mat_df$b, pb_mat_df$b, pb_mat_df$b)
  
  if (nrow(pb_mat_df) > 0) {
    for (i in seq_len(nrow(pb_mat_df))) {
      p <- p %>% add_annotations(
        x = pb_mat_df$p[i],
        y = pb_mat_df$b[i],
        xshift = layouts_xshift[[1]][i],
        yshift = layouts_yshift[[1]][i],
        text = pb_mat_df$self_constr[i],
        hoverinfo = "skip",
        font = list(size = 12 * text_size, color = "black"),
        showarrow = FALSE,
        xanchor = layouts_xanchor[[1]][i],
        yanchor = layouts_yanchor[[1]][i]
      )
    }
  }

  # plotly_build() resolves each add_annotations() call above into TWO
  # identical entries in the built layout (same text, same position) -
  # a quirk of chaining add_annotations() in a loop after add_markers() in
  # this plotly version. Building the plot ourselves here and keeping only
  # the first copy of each label avoids shipping a widget whose annotation
  # list is silently double the size it should be.
  p <- plotly::plotly_build(p)
  ann <- p$x$layout$annotations
  if (length(ann) > 0) {
    seen_ann_text <- character(0)
    keep <- vapply(ann, function(a) {
      if (a$text %in% seen_ann_text) return(FALSE)
      seen_ann_text[[length(seen_ann_text) + 1]] <<- a$text
      TRUE
    }, logical(1))
    p$x$layout$annotations <- ann[keep]
  }

  # Add Javascript interactivity
  pb_data <- list(
    dict = t,
    constructs = pb_mat_df$self_constr,
    col_rg = .construct_colors(wimp, mode = "red/green")[, "color"],
    col_gs = .construct_colors(wimp, mode = "grey scale")[, "color"],
    col_cb = .construct_colors(wimp, mode = "colorblind")[, "color"],
    col_dk = .construct_colors(wimp, mode = "dark")[, "color"],
    col_pt = .construct_colors(wimp, mode = "pastel")[, "color"],
    col_vd = .construct_colors(wimp, mode = "viridis")[, "color"],
    text_size = text_size,
    orig_x = pb_mat_df$p,
    orig_y = pb_mat_df$b,
    layouts_x = layouts_x,
    layouts_y = layouts_y,
    layouts_xshift = layouts_xshift,
    layouts_yshift = layouts_yshift,
    layouts_xanchor = layouts_xanchor,
    layouts_yanchor = layouts_yanchor
  )

  js_pb_panel <- "
    function(el, p_x, data) {
      var x = data;
      var pbUserAdjustedTextSize = false;
      var settingsModal = document.createElement('div');
      settingsModal.id = 'pb_settings_modal';
      Object.assign(settingsModal.style, {
        position: 'absolute', top: '10px', right: '10px',
        width: '90%', maxWidth: '300px', maxHeight: '80vh', overflowY: 'auto', boxSizing: 'border-box',
        backgroundColor: '#fff', zIndex: '2000', padding: '20px', borderRadius: '8px', 
        boxShadow: '0 4px 20px rgba(0,0,0,0.2)', border: '1px solid #eaeaea', display: 'none', 
        fontFamily: 'Inter, Roboto, sans-serif'
      });
      
      var pbHTML = '<div style=\"display:flex; justify-content:space-between; align-items:center; border-bottom:1px solid #eaeaea; padding-bottom:10px; margin-bottom:15px;\">' +
                   '<h3 style=\"margin:0; color:#444; font-size:14px;\">' + (x.dict.vis_options || 'Ajustes') + '</h3>' +
                   '<span id=\"close_pb_settings\" style=\"cursor:pointer; font-size:20px; font-weight:bold; color:#888; line-height:1;\">&times;</span>' +
                   '</div>';
                   
      pbHTML += '<div style=\"margin-bottom:15px;\">' +
                '<label style=\"display:block; margin-bottom:5px; font-weight:bold; color:#444; font-size:13px;\">' + (x.dict.color_palette || 'Paleta') + '</label>' +
                '<select id=\"pb_palette_sel\" style=\"width:100%; padding:4px; border-radius:4px;\">' +
                '<option value=\"rg\" selected>' + (x.dict.pal_redgreen || 'Red-Green') + '</option>' +
                '<option value=\"cb\">' + (x.dict.pal_colorblind || 'Colorblind') + '</option>' +
                '<option value=\"gs\">' + (x.dict.pal_greyscale || 'Greyscale') + '</option>' +
                '<option value=\"dk\">' + (x.dict.pal_dark || 'Dark') + '</option>' +
                '<option value=\"col_pt\">' + (x.dict.pastel || 'Pastel') + '</option>' +
                '<option value=\"col_vd\">' + (x.dict.viridis || 'Viridis') + '</option>' +
                '</select></div>';
                
      pbHTML += '<div style=\"margin-bottom:15px;\">' +
                '<label style=\"display:block; margin-bottom:5px; font-weight:bold; color:#444; font-size:13px;\">Tamaño del Texto</label>' +
                '<input type=\"range\" id=\"pb_text_size\" min=\"0.5\" max=\"2.5\" step=\"0.1\" value=\"' + x.text_size + '\" style=\"width:100%; accent-color:#8cc63f;\">' +
                '</div>';
                
      pbHTML += '<div style=\"margin-bottom:15px;\">' +
                '<label style=\"display:block; margin-bottom:5px; font-weight:bold; color:#444; font-size:13px;\">' + (x.dict.filter_constructs || 'Filtrar Constructos') + '</label>' +
                '<div id=\"pb_filter_list\" style=\"max-height:180px; overflow-y:auto; border:1px solid #ddd; padding:5px; border-radius:4px; font-size:12px; background:#f9f9f9;\"></div>' +
                '</div>';
                
      pbHTML += '<div style=\"margin-bottom:15px;\">' +
                '<label style=\"display:block; margin-bottom:5px; font-weight:bold; color:#444; font-size:13px;\">Etiquetas de Constructos</label>' +
                '<div style=\"margin-bottom:5px; display:flex; gap:10px;\">' +
                '<button id=\"btn_shuffle_labels\" style=\"flex:1; padding:6px; background:#f0f0f0; border:1px solid #ccc; border-radius:4px; cursor:pointer; font-size:12px; transition:0.2s;\">Reordenar</button>' +
                '<button id=\"btn_manual_adj\" style=\"flex:1; padding:6px; background:#f0f0f0; border:1px solid #ccc; border-radius:4px; cursor:pointer; font-size:12px; transition:0.2s;\">Ajuste Manual</button>' +
                '</div></div>';
                
      settingsModal.innerHTML = pbHTML;
      var container = el.closest('.wct-tab-content') || el.parentElement;
      container.appendChild(settingsModal);
      
      var filterContainer = settingsModal.querySelector('#pb_filter_list');
      x.constructs.forEach(function(lbl, idx) {
         var div = document.createElement('div');
         div.style.marginBottom = '4px';
         div.innerHTML = '<label style=\"cursor:pointer; display:flex; align-items:center; color:#555;\"><input type=\"checkbox\" checked value=\"' + idx + '\" class=\"pb-construct-cb\" style=\"margin-right:6px; accent-color:#8cc63f;\"> ' + lbl + '</label>';
         filterContainer.appendChild(div);
      });
      
      settingsModal.querySelector('#close_pb_settings').onclick = function() { settingsModal.style.display = 'none'; };
      
      // Inject Settings Button into the existing flexbox
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
         if (x.layouts_x && x.layouts_x.length > 0) {
             currentLayout = (currentLayout + 1) % x.layouts_x.length;
             updatePB();
         }
      };
      
      var updatePB = function() {
        var pal = settingsModal.querySelector('#pb_palette_sel').value;
        var txtSz = parseFloat(settingsModal.querySelector('#pb_text_size').value);
        
        var label_txt_sz = 'Text Size';
        if (x.dict.text_size) label_txt_sz = x.dict.text_size;
        else if (x.dict.hm_settings && x.dict.hm_settings.includes('Ajustes')) label_txt_sz = 'Tamaño del Texto';
        settingsModal.querySelector('#pb_text_size').previousSibling.textContent = label_txt_sz;
        
        var activeIndices = [];
        settingsModal.querySelectorAll('.pb-construct-cb').forEach(function(cb) {
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
           new_t.push('Construct: ' + x.constructs[idx] + '<br>P: ' + x.orig_x[idx] + '<br>B: ' + x.orig_y[idx]);
        }
        
        var restyleData = {
           x: [new_x],
           y: [new_y],
           'marker.color': [new_c],
           text: [new_t]
        };
        Plotly.restyle(el, restyleData, [0]);
        
        // Position the surviving labels with the SAME real-size-aware
        // layout fitPbLabels() uses (computePbLayout, defined below), not
        // the original R-computed positions (origAnnotations/x.layouts_*) -
        // those were solved for an 800x600 reference canvas and for every
        // construct, not the plot's actual current pixel size or whatever
        // subset the filter checkboxes leave active. Falling back to them
        // here is exactly what made labels visibly jump back whenever the
        // text-size slider or a filter checkbox was touched.
        var seenTexts = {};
        var uniqueAnn = origAnnotations.filter(function(a) {
          if (!a.text || seenTexts[a.text]) return false;
          seenTexts[a.text] = true;
          return true;
        });
        var activeSet = {};
        activeIndices.forEach(function(idx) { activeSet[x.constructs[idx]] = true; });
        var activeAnn = uniqueAnn.filter(function(a) { return activeSet[a.text]; });

        var fontPx = Math.max(6, Math.round(12 * txtSz));
        var w = el.clientWidth || 640;
        var marginL = el._fullLayout ? el._fullLayout.margin.l : 55;
        var marginB = el._fullLayout ? el._fullLayout.margin.b : 40;
        var innerW = Math.max(20, w - marginL - 20);
        var innerH = Math.max(20, (el.clientHeight || 480) - marginB - 10);
        var activeX = activeAnn.map(function(a) {
          var k = x.constructs.indexOf(a.text); return x.orig_x[k];
        });
        var activeY = activeAnn.map(function(a) {
          var k = x.constructs.indexOf(a.text); return x.orig_y[k];
        });
        var fresh = activeAnn.length
          ? computePbLayout(activeX, activeY, activeAnn.map(function(a) { return a.text; }), innerW, innerH, fontPx)
          : {xshift: [], yshift: [], xanchor: [], yanchor: []};

        var newAnnotations = activeAnn.map(function(a, i) {
          var ann = JSON.parse(JSON.stringify(a));
          ann.font.size = fontPx;
          if (fresh.xshift[i] !== undefined) {
            ann.xshift = fresh.xshift[i]; ann.yshift = fresh.yshift[i];
            ann.xanchor = fresh.xanchor[i]; ann.yanchor = fresh.yanchor[i];
          }
          return ann;
        });
        Plotly.relayout(el, { annotations: newAnnotations });
      };
      
      settingsModal.querySelector('#pb_palette_sel').onchange = updatePB;
      settingsModal.querySelector('#pb_text_size').oninput = function() { pbUserAdjustedTextSize = true; updatePB(); };
      settingsModal.querySelectorAll('.pb-construct-cb').forEach(function(cb) {
         cb.onchange = updatePB;
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

      // The collision-free placement further up (.calculate_pb_layouts, in
      // R) solves overlap in a FIXED 800x600 reference canvas - correct
      // only when the real plot happens to render near that size. Once the
      // widget is resized to anything else, the real pixel spacing between
      // points no longer matches what that layout planned around: two
      // points comfortably apart in an 800px-wide reference end up much
      // closer together in a real 300px-wide plot, so their (roughly
      // fixed-size) label boxes collide even though the precomputed
      // xanchor/xshift were technically non-overlapping on paper. The fix
      // is to re-run placement here in JS, in the plot's ACTUAL current
      // pixel dimensions, every time it resizes - not just rescale the font
      // on positions computed for a different size.
      //
      // Two earlier attempts here: a spring that let labels drift far from
      // their point (looked broken), then one constrained to spin at a
      // near-fixed radius (too busy with this many labels). Back to 8
      // fixed slots around the point - right/left/top/bottom and the four
      // diagonals, same small offsets a human would pick by hand - but
      // instead of just counting overlapping pixels to score each slot,
      // score it the way same-sign charges would feel it: every other
      // point and every already-placed label contributes a repulsion
      // energy that grows sharply at short range (inverse-square), so a
      // slot that's merely close to a crowd loses out to one with real
      // breathing room, not only one that's technically overlap-free.
      function computePbLayout(xs, ys, labels, canvasW, canvasH, fontPx) {
        var n = xs.length;
        var xMin = Math.min.apply(null, xs), xMax = Math.max.apply(null, xs);
        var yMin = Math.min.apply(null, ys), yMax = Math.max.apply(null, ys);
        var xRange = (xMax - xMin) || 1, yRange = (yMax - yMin) || 1;
        // +y here means up, matching how Plotly's own xshift/yshift work
        // (a positive yshift moves an annotation up on screen) - so these
        // can be used as shift values directly, with no sign flip.
        var px = xs.map(function(v) { return (v - xMin) / xRange * canvasW; });
        var py = ys.map(function(v) { return (v - yMin) / yRange * canvasH; });
        var charW = fontPx * 0.58, charH = fontPx * 1.25;
        var labelW = labels.map(function(l) { return (l ? String(l).length : 4) * charW; });
        var labelH = labelW.map(function() { return charH; });

        var gapOrtho = 5, gapDiag = 2;
        var candidates = [
          {dx: gapOrtho, dy: 0, xa: 'left', ya: 'middle'},
          {dx: -gapOrtho, dy: 0, xa: 'right', ya: 'middle'},
          {dx: 0, dy: gapOrtho, xa: 'center', ya: 'bottom'},
          {dx: 0, dy: -gapOrtho, xa: 'center', ya: 'top'},
          {dx: gapDiag, dy: gapDiag, xa: 'left', ya: 'bottom'},
          {dx: -gapDiag, dy: gapDiag, xa: 'right', ya: 'bottom'},
          {dx: gapDiag, dy: -gapDiag, xa: 'left', ya: 'top'},
          {dx: -gapDiag, dy: -gapDiag, xa: 'right', ya: 'top'}
        ];

        var placed = []; // {cx, cy, left, right, bottom, top} of labels placed so far
        var out = {xshift: [], yshift: [], xanchor: [], yanchor: []};

        for (var i = 0; i < n; i++) {
          var tw = labelW[i], th = labelH[i];
          var bestCand = candidates[0], bestEnergy = Infinity, bestOverlap = Infinity;

          for (var c = 0; c < candidates.length; c++) {
            var cand = candidates[c];
            var cx = px[i] + cand.dx, cy = py[i] + cand.dy;
            var left = cx, right = cx;
            if (cand.xa === 'left') { left = cx; right = cx + tw; }
            else if (cand.xa === 'right') { left = cx - tw; right = cx; }
            else { left = cx - tw / 2; right = cx + tw / 2; }
            var bottom, top;
            if (cand.ya === 'bottom') { bottom = cy; top = cy + th; }
            else if (cand.ya === 'top') { bottom = cy - th; top = cy; }
            else { bottom = cy - th / 2; top = cy + th / 2; }
            var boxCx = (left + right) / 2, boxCy = (bottom + top) / 2;

            // Hard overlap count still breaks ties first - a slot that
            // overlaps nothing beats one with lower energy but a visible
            // collision.
            var overlap = 0;
            for (var j = 0; j < n; j++) {
              if (j === i) continue;
              if (px[j] >= left && px[j] <= right && py[j] >= bottom && py[j] <= top) overlap += 1;
            }
            for (var b = 0; b < placed.length; b++) {
              var box = placed[b];
              if (!(left >= box.right || right <= box.left || bottom >= box.top || top <= box.bottom)) overlap += 1;
            }

            // Electrostatic-style potential energy at this slot's box
            // center: every other point and every placed label pushes back
            // harder the closer they are.
            var energy = 0;
            for (var k = 0; k < n; k++) {
              if (k === i) continue;
              var dk = Math.max(6, Math.hypot(boxCx - px[k], boxCy - py[k]));
              energy += 1 / (dk * dk);
            }
            for (var b2 = 0; b2 < placed.length; b2++) {
              var db = Math.max(6, Math.hypot(boxCx - placed[b2].cx, boxCy - placed[b2].cy));
              energy += 30 / (db * db); // another label crowds a slot more than a bare point does
            }

            if (overlap < bestOverlap || (overlap === bestOverlap && energy < bestEnergy)) {
              bestOverlap = overlap; bestEnergy = energy; bestCand = cand;
            }
          }

          out.xshift.push(bestCand.dx); out.yshift.push(bestCand.dy);
          out.xanchor.push(bestCand.xa); out.yanchor.push(bestCand.ya);

          var fcx = px[i] + bestCand.dx, fleft, fright, fbottom, ftop;
          if (bestCand.xa === 'left') { fleft = fcx; fright = fcx + tw; }
          else if (bestCand.xa === 'right') { fleft = fcx - tw; fright = fcx; }
          else { fleft = fcx - tw / 2; fright = fcx + tw / 2; }
          var fcy = py[i] + bestCand.dy;
          if (bestCand.ya === 'bottom') { fbottom = fcy; ftop = fcy + th; }
          else if (bestCand.ya === 'top') { fbottom = fcy - th; ftop = fcy; }
          else { fbottom = fcy - th / 2; ftop = fcy + th / 2; }
          placed.push({cx: (fleft + fright) / 2, cy: (fbottom + ftop) / 2, left: fleft, right: fright, bottom: fbottom, top: ftop});
        }
        return out;
      }

      // Auto-shrink the construct labels, dots and axis titles as the
      // widget gets narrower, and recompute their placement for the real
      // current size instead of just rescaling font on a layout planned
      // for a different one. Skipped once the user has manually set a text
      // size, so it never fights their choice.
      var pbBaseAnnotations = JSON.parse(JSON.stringify(el.layout.annotations || []));
      // The construct labels come back doubled here (same text, same
      // position, twice) - keep only the first copy of each. Without this,
      // our resize logic only had 21 fresh positions (one per construct)
      // for 42 annotations, so the second copy of each label never got
      // updated and was left behind at its original spot: what looked like
      // a duplicated label was really its stale twin staying put while
      // the real one moved.
      (function() {
        var seen = {};
        pbBaseAnnotations = pbBaseAnnotations.filter(function(a) {
          if (seen[a.text]) return false;
          seen[a.text] = true;
          return true;
        });
      })();
      function fitPbLabels() {
        if (pbUserAdjustedTextSize) return;
        var w = el.clientWidth; if (!w) return;
        var scale = Math.max(0.55, Math.min(1, w / 640));
        var fontPx = Math.max(8, Math.round(12 * scale));
        var marginL = Math.max(40, Math.round(55 * scale));
        var marginB = Math.max(28, Math.round(40 * scale));
        var innerW = Math.max(20, w - marginL - 20);
        var innerH = Math.max(20, (el.clientHeight || 480) - marginB - 10);
        var fresh = computePbLayout(x.orig_x, x.orig_y, x.constructs, innerW, innerH, fontPx);
        var anns = pbBaseAnnotations.map(function(a, i) {
          var b = Object.assign({}, a);
          b.font = Object.assign({}, a.font, {size: fontPx});
          if (fresh.xshift[i] !== undefined) {
            b.xshift = fresh.xshift[i]; b.yshift = fresh.yshift[i];
            b.xanchor = fresh.xanchor[i]; b.yanchor = fresh.yanchor[i];
          }
          return b;
        });
        var titleFontPx = Math.max(11, Math.round(20 * scale));
        var standoffPx = Math.max(2, Math.round(4 * scale));
        var tickFontPx = Math.max(9, Math.round(13 * scale));
        Plotly.relayout(el, {
          annotations: anns,
          'xaxis.title.font.size': titleFontPx, 'xaxis.title.standoff': standoffPx,
          'yaxis.title.font.size': titleFontPx, 'yaxis.title.standoff': standoffPx,
          'xaxis.tickfont.size': tickFontPx, 'yaxis.tickfont.size': tickFontPx,
          'margin.l': marginL, 'margin.b': marginB
        }).then(clampPbEdgeLabels);
        Plotly.restyle(el, {'marker.size': Math.max(4, Math.round(7 * scale))}, [0]);
      }
      // A label anchored to extend leftward from a point close to the left
      // edge (x = 0, where Presence is defined to start) - or rightward
      // from a point near the right edge - can render past that edge and
      // get clipped, especially at narrow widths where the precomputed
      // collision-free layout has less room to work with than it assumed.
      // Once real positions are on screen, re-anchor any label whose
      // rendered edge crosses the plot's actual drawing area so it points
      // back inward instead.
      function clampPbEdgeLabels() {
        var elRect = el.getBoundingClientRect();
        var plotLeftPx = elRect.left + el._fullLayout.margin.l;
        var plotRightPx = elRect.right - el._fullLayout.margin.r;
        var xrange = el._fullLayout.xaxis.range;
        var xSpan = xrange[1] - xrange[0];
        var byText = {};
        el.querySelectorAll('.annotation-text').forEach(function(n) { byText[n.textContent] = n; });
        var changed = false;
        var anns = (el.layout.annotations || []).map(function(a) {
          var b = Object.assign({}, a);
          var node = byText[a.text];
          if (!node) return b;
          var r = node.getBoundingClientRect();
          // Only flip a label toward whichever edge its OWN point is
          // actually near - re-anchoring a label that overflows simply
          // because it's long, with its point nowhere near that edge,
          // would just swing it past the opposite edge instead.
          var frac = (a.x - xrange[0]) / xSpan;
          if (r.left < plotLeftPx - 1 && b.xanchor !== 'left' && frac < 0.5) {
            b.xanchor = 'left'; b.xshift = 3; changed = true;
          } else if (r.right > plotRightPx + 1 && b.xanchor !== 'right' && frac > 0.5) {
            b.xanchor = 'right'; b.xshift = -3; changed = true;
          }
          return b;
        });
        if (changed) Plotly.relayout(el, {annotations: anns});
      }
      var pbLastWidth = el.clientWidth;
      setTimeout(fitPbLabels, 50);
      setInterval(function() {
        var w = el.clientWidth;
        if (w && w !== pbLastWidth) { pbLastWidth = w; fitPbLabels(); }
      }, 300);

      // Plotly's tickformat can't drop a leading zero on its own (\".2\"
      // instead of \"0.2\"), so strip it straight from the rendered tick
      // text - this fires on every redraw, including the relayout calls
      // above, so it stays correct as the plot rescales.
      function stripPbLeadingZeros() {
        el.querySelectorAll('.xtick text, .ytick text').forEach(function(t) {
          var s = t.textContent;
          // Plotly renders the minus sign as U+2212, not a plain hyphen.
          var m = /^([-\\u2212]?)0(\\.\\d+)$/.exec(s);
          if (!m) return;
          var isZero = /^\\.0+$/.test(m[2]);
          var out = isZero ? '0' : (m[1] + m[2]);
          if (out !== s) t.textContent = out;
        });
      }
      el.on('plotly_afterplot', stripPbLeadingZeros);
      setTimeout(stripPbLeadingZeros, 60);
    }
  "
  
  p <- p %>% htmlwidgets::onRender(js_pb_panel, data = pb_data)

  return(p)
}
