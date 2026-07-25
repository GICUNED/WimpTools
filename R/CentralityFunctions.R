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

  # Normalize coordinates for the label optimization algorithm
  # Map from data coordinates to normalized [0,1] range
  x_range <- range(pb_mat_df$p)
  y_range <- range(pb_mat_df$b)
  x_span <- diff(x_range)
  y_span <- diff(y_range)

  # Add small padding to avoid edge issues
  x_padding <- x_span * 0.1
  y_padding <- y_span * 0.1

  norm_x <- (pb_mat_df$p - (x_range[1] - x_padding)) / (x_span + 2 * x_padding)
  norm_y <- (pb_mat_df$b - (y_range[1] - y_padding)) / (y_span + 2 * y_padding)

  # Apply smart label positioning algorithm
  optimized_positions <- .smart_label_positions(
    x_coords = norm_x,
    y_coords = norm_y,
    labels = pb_mat_df$self_constr,
    distance = 8,
    text_size = 12 * text_size
  )

  # Convert optimized positions back to original data scale
  if (nrow(optimized_positions) > 0) {
    # Convert normalized coordinates back to data coordinates
    optimized_positions$x_data <- optimized_positions$x *
      (x_span + 2 * x_padding) + (x_range[1] - x_padding)
    optimized_positions$y_data <- optimized_positions$y *
      (y_span + 2 * y_padding) + (y_range[1] - y_padding)
  }

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
      marker = list(color = pb_mat_df$color, size = 7,
                    line = list(color = "black", width = 1)),
      text = ~paste("Construct:", self_constr, "<br>P:", p, "<br>B:", b),
      hoverinfo = "text"
    )

  # Add optimized annotations
  if (nrow(optimized_positions) > 0) {
    for (i in seq_len(nrow(optimized_positions))) {
      p <- p %>% add_annotations(
        x = pb_mat_df$p[i],
        y = pb_mat_df$b[i],
        text = optimized_positions$label[i],
        hoverinfo = "skip",
        font = list(size = 12 * text_size, color = "black"),
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

  return(p)
}
