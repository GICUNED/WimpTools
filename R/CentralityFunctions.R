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
                                standoff = 25)),
      yaxis = list(title = list(text = t$pb_axis_y,
                                font = list(size = 20), standoff = 25)),
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
        
        var newAnnotations = [];
        var seenTexts = {};
        for (var i = 0; i < origAnnotations.length; i++) {
           var ann = JSON.parse(JSON.stringify(origAnnotations[i]));
           if (!ann.text || seenTexts[ann.text]) continue;
           seenTexts[ann.text] = true;
           
           // Find the construct index by matching text
           var c_idx = -1;
           for (var k = 0; k < x.constructs.length; k++) {
               if (x.constructs[k] === ann.text) {
                   c_idx = k;
                   break;
               }
           }
           
           if (c_idx !== -1 && x.layouts_x && x.layouts_x.length > currentLayout) {
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
               ann.font.size = 12 * txtSz;
               newAnnotations.push(ann);
           }
        }
        Plotly.relayout(el, { annotations: newAnnotations });
      };
      
      settingsModal.querySelector('#pb_palette_sel').onchange = updatePB;
      settingsModal.querySelector('#pb_text_size').oninput = updatePB;
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
    }
  "
  
  p <- p %>% htmlwidgets::onRender(js_pb_panel, data = pb_data)

  return(p)
}
