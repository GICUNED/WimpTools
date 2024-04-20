## CENTRALITY FUNCTIONS ##

# PH Indices ------------------------------------------------------------

#' Calculate Presence and Hierarchy Indices
#'
#' This function computes the presence (P, frequency of occurrence) and
#' hierarchy (H, influence on others) indices for constructs within an implication grid.
#' It can standardize these indices based on the maximum degree if required.
#'
#' @param wimp An object of class 'wimp', which contains an implication grid
#'   and associated constructs.
#' @param method A character string specifying the method used to calculate
#'   the degree indices. Default is "weight". Other methods may be available
#'   depending on the implementation of 'degree_index' function.
#' @param std A logical value indicating whether to standardize the P and H
#'   indices based on the maximum total degree of the constructs. Defaults to FALSE
#'
#' @return A 2-column matrix with the presence index 'p' and the hierarchy index 'h'
#'   for each construct. If 'std' is TRUE, these values are standardized.
#'
#' @export
#'
#' @examples
#' # Assuming 'wimp' is an object with an implication grid
#' ph_indices <- ph_index(wimp)
#' ph_indices_std <- ph_index(wimp, std = TRUE)
#' ph_indices_wnorm_non_std <- ph_index(wimp, method = "wnorm", std = FALSE)
#'

ph_index <- function(wimp, method = "wnorm", std = 'none'){

  # P-H calculation--------------
  # Connectivity of constructs
  c.io <- degree_index(wimp, method = method)

  # Rearrange In - Out columns
  c.io <- c.io[, c(2,1,3)]

  # Extract In and Out columns
  in.out <- c.io[, 1:2]

  # Linear Transformation matrix
  coef <- 1 / sqrt(2)
  coef.matrix <- matrix(c(coef, -coef, coef, coef), nrow = 2)

  # Calculate P - H matrix
  ph.mat <- in.out %*% t(coef.matrix)
  colnames(ph.mat) <- c("p", "h")

  # Standardization--------------
  if (std == 'vertices'){
    vertices <- length(wimp$constructs$constructs)

    coef.max.p <- 2*coef*(vertices-1)
    coef.max.h <- coef*(vertices-1)

    ph.mat[,"p"] <- ph.mat[,"p"]/coef.max.p
    ph.mat[,"h"] <- ph.mat[,"h"]/coef.max.h

  }else if (std == 'edges'){
    c.direct.io <- degree_index(wimp, method = "simple")
    c.direct.io <- c.direct.io[, c(2,1,3)]

    edges <- sum(c.direct.io[, 2]) # The number of edges is the sum of all outgoing connections

    coef.max.p <- edges * coef
    coef.max.h <- edges * coef

    ph.mat[,"p"] <- ph.mat[,"p"]/coef.max.p
    ph.mat[,"h"] <- ph.mat[,"h"]/coef.max.h

  }else if (std == 'max_edges'){
    edges <- max((c.io[, 2])) # Max outgoing connections

    coef.max.p <- edges
    coef.max.h <- edges

    ph.mat[,"p"] <- ph.mat[,"p"]/coef.max.p
    ph.mat[,"h"] <- ph.mat[,"h"]/coef.max.h


  }else if (std == "density"){
    vertices <- length(wimp$constructs$constructs)

    max.edges <- vertices * (vertices -1) # Maximum theoretical number of edges

    c.direct.io <- degree_index(wimp, method = "simple")
    c.direct.io <- c.direct.io[, c(2,1,3)]
    total.edges <- sum(c.direct.io[, 2]) # The number of edges is the sum of all outgoing connections

    dens <- total.edges/ max.edges

    ph.mat[,"p"] <- ph.mat[,"p"]*dens
    ph.mat[,"h"] <- ph.mat[,"h"]*dens

  }

  return(ph.mat)
}


# Mahalanobis Indices ------------------------------------------------------------

#' Calculate Mahalanobis Indices for Constructs
#'
#' This function calculates the Mahalanobis distance for each construct
#' in a given PH matrix obtained from a `wimp` object. It also determines
#' whether each construct is considered "central" based on a chi-square
#' cutoff, which is calculated using a specified significance level.
#'
#' @param wimp A `wimp` object containing the data from which the PH matrix
#'   is derived.
#' @param method A character string specifying the method used for calculating
#'   the PH matrix. Default is `"wnorm"`.
#' @param std A logical value indicating whether the data should be
#'   standardized before calculating the Mahalanobis distance. Default is `FALSE`.
#'
#' @return A matrix that includes the original PH matrix with two additional
#'   columns: one for the Mahalanobis distance of each construct and another
#'   indicating whether the construct is considered "central" based on the
#'   chi-square cutoff.
#'
#' @export
#'
#' @examples
#' result <- mahalanobis_index(wimp, method = "weight", std = FALSE)
#' head(result)
#'

mahalanobis_index <- function(wimp, method = "wnorm", std = 'none', sign.level = 0.2){

  # Obtain PH matrix associated to the wimp object
  ph.mat <- ph_index(wimp = wimp, method = method, std = std)

  # Greater weight is given to a positive hierarchy value
  #ph.mat[,"h"] <- ifelse(ph.mat[,"h"] >0, ph.mat[,"h"]^2, ph.mat[,"h"])

  # Mahalanobis distance--------

  # Covariance matrix
  cov.matrix <- cov(ph.mat)
  # Means vector
  means.vect <- colMeans(ph.mat)
  # Calculate the Mahalanobis distance for each observation in the dataframe
  ph.mat.df <- as.data.frame(ph.mat)
  m.dist <- mahalanobis(ph.mat.df, center = means.vect, cov = cov.matrix)

  # Add column to the results matrix
  phm.mat <- cbind(ph.mat, m.dist)

  # Determine the chi-square cutoff for the given significance level--------------
  # Degrees of freedom
  df <- ncol(ph.mat)
  # Chi-Square Cutoff
  chi.square.cutoff <- qchisq(1 - sign.level, df)

  # Annotate observations as "hub" constructs based on the chi-square cutoff
  # and being "non-superficial" constructs on the P axis
  #----------------------
  # Mean and standard deviation of the distribution
  mean.p <- mean(ph.mat.df$p)
  sd.p <- sd(ph.mat.df$p)

  # Cutoff point - Defined with respect to distribution
  p.cut <- qnorm(0.15, mean = mean.p, sd = sd.p)

  # Filter out P values that are less than the cutoff point
  hub <- ifelse(phm.mat[,"m.dist"] > chi.square.cutoff & phm.mat[,"p"] > p.cut, TRUE, FALSE)
  #----------------------

  # Add column to the results matrix
  phmc.mat <- cbind(phm.mat, hub)
  return(phmc.mat)
}

# PCA Indices ------------------------------------------------------------

#' Calculate centrality scores based on Principal Component Analysis (PCA)
#'
#' This function computes centrality scores for constructs in a WIMP object using
#' PCA on a specified matrix (direct, weights, or implications). It returns the centrality scores
#' based on the loadings of the specified number of principal components.
#'
#' @param wimp A list object representing the WIMP structure, which must include
#'   a scores list with matrices: direct, weights, and implications.
#' @param matrix A character string specifying which wimp matrix to use for PCA.
#'   Valid options are "direct", "weights", or "implications".
#' @param pr.comp Integer indicating the number of principal components to use
#'   for calculating centrality scores. Defaults to 2.
#'
#' @return A dataframe with construct names and their centrality scores.
#'
#' @examples
#' pca_index(my_wimp, matrix = "weights", pr.comp = 2)
#'
#' @export
#'
#' @throws Error if method is not one of the allowed values. Also throws an error if pr.comp
#'   exceeds the available number of principal components.
#'

pca_index <- function(wimp, matrix = "implications", pr.comp = 2){
  # Validate the specified method and its presence in wimp$scores
  if (!matrix %in% c("direct", "weights", "implications")) {
    stop("El método especificado debe ser 'direct', 'weights' o 'implications'.")
  }

  # Access the matrix based on the specified method
  adj.matrix <- wimp$scores[[matrix]]

  pca.result <- NULL
  # Attempt to perform PCA. Catch potential calc exceptions
  pca.result <- tryCatch({
    prcomp(adj.matrix, center = TRUE, scale = TRUE)
  }, warning = function(w) {
    message("Advertencia: ", w$message)
  }, error = function(e) {
    message("Error en el cálculo del PCA: ", e$message)
  })

  # Check if PCA was successful
  if (is.null(pca.result)) {
    # PCA failed, assign NA to all centrality values
    centrality <- rep(NA, length(wimp$constructs$constructs))
  } else {
    # PCA succeeded, calculate variance explained and loadings
    explained.variance <- pca.result$sdev^2 / sum(pca.result$sdev^2)

    # Validate if pr.comp is within the allowable range
    if (pr.comp > length(explained.variance)) {
      stop("pr.comp excede el número de componentes principales disponibles.")
    }

    # Calculate loadings for the specified number of principal components
    loadings.pca <- pca.result$rotation
    loadings.add <- rowSums((loadings.pca[, 1:pr.comp]^2) * explained.variance[1:pr.comp])
    centrality <- abs(loadings.add)
  }

  # Create a dataframe with construct names and centrality scores from the sum of loadings in principal components
  pca.df <- data.frame(
    constructs = wimp$constructs$constructs,
    leftpoles = wimp$constructs$left.poles,
    rightpoles = wimp$constructs$right.poles,
    centrality = centrality
  )

  return(pca.df)
}


# Eigen Indices ------------------------------------------------------------

#' Calculate Centrality Using Eigenvalue Decomposition
#'
#' This function calculates centrality scores for constructs within a WIMP (Web-based Ideographic Measures Package) structure,
#' based on the eigenvalue decomposition of a specified adjacency matrix from the WIMP scores. It supports analyzing centrality
#' using the 'direct', 'weights', or 'implications' matrices.
#'
#' @param wimp wimp An object of class 'wimp' (weighted implications grid)
#' @param matrix A character string specifying which matrix to use for the centrality analysis. Accepted values are
#'        'direct', 'weights', or 'implications'. Default is 'implications'.
#' @param num.vectors An integer specifying the number of eigenvectors to use for computing centrality scores.
#' @return A dataframe containing the constructs' names and their respective centrality scores.
#'
#' @export
#'
#' @throws Error if method is not one of the allowed values. Also throws an error if num.vectors
#'   exceeds the available number of eigenvectors

eigen_index <- function(wimp, matrix = "implications", num.vectors = 2) {
  # Validate the specified method
  if (!matrix %in% c("direct", "weights", "implications")) {
    stop("method debe ser 'direct', 'weights' o 'implications'.")
  }

  # Access the matrix based on the specified method
  adj.matrix <- wimp$scores[[matrix]]

  # Compute eigenvectors and eigenvalues
  results <- eigen(adj.matrix)

  # Ensure the number of requested eigenvectors does not exceed available components
  if (num.vectors > length(results$values)) {
    stop("num.vectors excede el número de vectores propios disponibles.")
  }

  # Calculate centrality using the specified number of eigenvectors. This calculation takes the absolute value
  # of the real part of each eigenvector, squares it, and then multiplies it by the real part of the corresponding
  # eigenvalue to compute centrality scores.
  centralidad <- Reduce(`+`, lapply(1:num.vectors, function(i) {
    Re(results$vectors[, i])^2 * Re(results$values[i])
  }))

  # Create a dataframe with the centrality results
  df.centrality <- data.frame(
    constructs = wimp$constructs$constructs,
    leftpoles = wimp$constructs$left.poles,
    rightpoles = wimp$constructs$right.poles,
    centrality = abs(centralidad)
  )

  return(df.centrality)
}

# Graphically represent constructs ------------------------------------------------

#' Scatter Plot of Constructs in P-H Space
#'
#' This function generates a scatter plot of constructs in the P-H space,
#' where P represents Presence (frequency of the construct) and H represents Hierarchy
#' (influence of the construct). Hub constructs are highlighted in red color, and
#' peripheral constructs in another, facilitating their visual identification.
#'
#' @param phm.mat A matrix where each row represents a construct and contains
#'        the P and H coordinates of the construct, as well as an indication of whether the
#'        construct is hub (1) or not (0). The matrix must have row names,
#'        which are used to label the constructs in the graph.
#'
#' @param mark.nva Boolean value that specifies if non-viable areas are to be marked in the
#'        graphic. Non-viable areas are parts of the P-H space where constructs cannot logically exist.
#'
#' @param mark.hub Boolean value that specifies if hub constructs have to be highlighted.
#'        Hub constructs are depicted with a distinct color and symbol to differentiate them from
#'        peripheral constructs.
#'
#' @param show.points Boolean value that specifies whether points should be displayed or not
#'
#' @return A Plotly object representing the generated scatter plot.
#'
#' @export
#'
#' @examples
#' graph_ph(phm.mat, mark.nva = TRUE, mark.hub = TRUE)

graph_ph <- function(..., mark.nva = TRUE, mark.hub = TRUE, show.points = TRUE) {

  # Extract the Mahalobis distance matrix for the given wimp
  phm.mat <- mahalanobis_index(...)
  # Convert the matrix to a dataframe
  phm.mat.df <- as.data.frame(phm.mat)
  # Assign the names of constructs from the row names of the matrix
  phm.mat.df$constructo <- rownames(phm.mat)
  # Assign the names of constructs in "Self"
  phm.mat.df$self.constr <- wimp$constructs$self.poles
  # Limits for the graph by the largest value of P or H dimensions. We add a small margin
  limit <- max(abs(phm.mat.df$p), abs(phm.mat.df$h)) * 1.1

  # Round P and H values
  phm.mat.df$p <- round(phm.mat.df$p, 3)
  phm.mat.df$h <- round(phm.mat.df$h, 3)

  # Add a new column for label color
  if (mark.hub)
    phm.mat.df$label.color <- ifelse(phm.mat.df$hub == 1, 'red', 'black')
  else
    phm.mat.df$label.color <- 'black'

  # Shapes for non-viable area if requested
  shapes <- if (mark.nva) {
    list(
      list(type = "path", path = paste("M 0,0 L", limit, ",", limit, " L0,", limit, " Z"),
           fillcolor = "palegreen", opacity = 0.2, line = list(color = "palegreen")),
      list(type = "path", path = paste("M 0,0 L", limit, ",", -limit, " L0,", -limit, " Z"),
           fillcolor = "palegreen", opacity = 0.2, line = list(color = "palegreen")),
      list(type = "line", x0 = 0, y0 = 0, x1 = limit, y1 = limit,
           xref = "x", yref = "y", line = list(color = "darkgreen", width = 1, dash = "dash")),
      list(type = "line", x0 = 0, y0 = 0, x1 = limit, y1 = -limit,
           xref = "x", yref = "y", line = list(color = "darkgreen", width = 1, dash = "dash"))
    )
  } else {
    list()  # No shapes if mark.nva is FALSE
  }

  # Initialize Plotly graph
  p <- plot_ly()

  # Set the layout of the graph
  p <- p %>%
    layout(title = 'Gráfica de Dispersión de Constructos en el Espacio P - H',
           xaxis = list(title = 'P - Presencialidad (frecuencia del constructo)'),
           yaxis = list(title = 'H - Jerarquía (influencia del constructo)'),
           plot_bgcolor = "white",
           font = list(family = "Arial"),
           showlegend = FALSE,
           shapes = shapes)

  # Add points or not based on show.points parameter
  if (show.points) {
    # Construct category colors
    phm.mat.df$color <- NA

    colors.mat <- construct_colors(wimp = wimp, mode = "red/green")
    phm.mat.df$color <- colors.mat[,"color"]

    if (mark.hub) {
      p <- p %>%
        add_markers(data = phm.mat.df, x = ~p, y = ~h,
                    marker = list(color = ~color, size = 10,
                                  symbol = ifelse(phm.mat.df$hub == 1, "star", "circle"),
                                  line = list(color = 'black', width = ifelse(phm.mat.df$hub == 1, 2, 1))),
                    text = ~paste('P:', p, '; H:', h), hoverinfo = 'text')
    } else {
      p <- p %>%
        add_markers(data = phm.mat.df, x = ~p, y = ~h,
                    #marker = list(color = colors, size = 10, line = list(color = 'black', width = 1)),
                    marker = list(color = phm.mat.df$color, size = 10, line = list(color = 'black', width = 1)),
                    text = ~paste('P:', p, '; H:', h), hoverinfo = 'text')
    }
  }

  # Add construct labels (annotations)
  if (mark.hub & !show.points) { # Labels are highlighted if centers are checked and no points are displayed
    p <- p %>%
      add_annotations(data = phm.mat.df[phm.mat.df$hub == 0, ], x = ~p, y = ~h, text = ~self.constr,
                      hovertext = ~paste('Constructo:', constructo, '\nP:', p, 'H:', h), hoverinfo = 'text',
                      font = list(size = 12, color = 'black'),
                      showarrow = FALSE, xanchor = 'center', yanchor = 'bottom',
                      yshift = 5)

    p <- p %>%
      add_annotations(data = phm.mat.df[phm.mat.df$hub == 1, ], x = ~p, y = ~h, text = ~self.constr,
                      hovertext = ~paste('Constructo:', constructo, '\nP:', p, 'H:', h), hoverinfo = 'text',
                      font = list(size = 12, color = 'darkgreen', family = "Arial Black, sans-serif",
                                  style = "normal"),
                      showarrow = FALSE, xanchor = 'center', yanchor = 'bottom',
                      yshift = 5)

  } else { # Labels are normal whether no centers are marked or points are displayed.
    p <- p %>%
      add_annotations(data = phm.mat.df, x = ~p, y = ~h, text = ~self.constr,
                      hovertext = ~paste('Constructo:', constructo, '\nP:', p, 'H:', h), hoverinfo = 'text',
                      font = list(size = 12, color = 'black'),
                      showarrow = FALSE, xanchor = 'center', yanchor = 'bottom',
                      yshift = 5)
  }

  return(p)
}


# Construct colors ------------------------------------------------------------

#' Construc color matrix based on its category
#'
#' This function computes the presence (P, frequency of occurrence) and
#' hierarchy (H, influence on others) indices for constructs within an implication grid.
#' It can standardize these indices based on the maximum degree if required.
#'
#' @param wimp An object of class 'wimp', which contains an implication grid
#'   and associated constructs.
#' @param mode Color mode for .color.selection (i.e, "red/green")
#'
#' @return A matrix containing construct colors based on its category and .color.selection.
#'
#' @export
#'
#' @examples
#'

construct_colors <- function(wimp, mode){
  # Construct category colors
  col.sel <- .color.selection(mode)
  color.mat <- matrix(data = 0, nrow = length(wimp$constructs$constructs), ncol = 1)
  rownames(color.mat) <- wimp$constructs$constructs
  colnames(color.mat) <- c("color")

  color.mat[wimp$constructs$discrepants,"color"] <- col.sel[1]
  color.mat[wimp$constructs$congruents,"color"] <- col.sel[2]
  color.mat[wimp$constructs$undefined,"color"] <- col.sel[3]
  color.mat[wimp$constructs$dilemmatic,"color"] <- col.sel[4]

  return(color.mat)
}

# Test optimal numbers of clusters ------------------------------------------------------------

#' Public test function for .optimal.num.clusters hide function
#'
#' This function computes the presence (P, frequency of occurrence) and
#'
#' @param
#'
#' @param
#'
#' @return
#'
#' @export
#'
#' @examples
#'

test_optimal_num_clusters <- function(...){
  # Invoke the hidden function
  return(.optimal.num.clusters(...))
}


# Dendrogram of constructs of a Wimp ------------------------------------------------------------
#' Constructs Dendrogram
#'
#' This function generates a dendrogram of constructs based on the Euclidean
#' distances in the P-H vector space. It automatically calculates the optimal
#' number of clusters for grouping the constructs.
#'
#' @param wimp A `wimp` object containing the constructs and the corresponding
#' scores in the P-H vector space. This object should have been generated using
#' the GridFCM.practicum package.
#'
#' @return A dendrogram visualizing the hierarchical clustering of constructs
#' based on their distances in the P-H space. The dendrogram highlights the
#' optimal clustering of constructs.
#'
#' @importFrom GridFCM.practicum ph_index
#'
#' @importFrom factoextra fviz_dend hcut
#'
#' @export
#'

constructs_dendrogram <- function(wimp){

  # Matrix of Euclidean distances between constructs in the vector space P-H
  ph.mat <- GridFCM.practicum::ph_index(wimp = wimp, method = "wnorm", std = FALSE)
  rownames(ph.mat) <- wimp$constructs$self.poles

  # Optimal number of clusters
  k<- .optimal.num.clusters(wimp)

  # Colors palette
  greens.palette <- c("#003300","#008000", "#3CB371")

  #Dendrogram
  dist.mat <- .mahalanobis.dist.matrix(ph.mat)
  #dist.mat <- dist(ph.mat, method = "euclidean")
  hclust.model <- hcut(dist.mat, k = k, method = "ward.D2", stand = TRUE, hc_func = "agnes")
  plot<- fviz_dend(hclust.model, rect = TRUE, cex = 0.7,
                   k_colors = greens.palette,horiz = TRUE,
                   main = 'Dendrograma de constructos por distancia en P-H')


  return(plot)
}

