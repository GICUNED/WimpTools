## CENTRALITY FUNCTIONS ##

# Degree Index Centrality -------------------------------------------------

#' Degree Index -- degree_index()
#'
#' @description Function to calculate the centrality of the constructs.
#'              In this case, centrality is understood as the degree of connection that each
#'              construct maintains with the rest, i.e. the number of links for each vertex.
#'
#' @param wimp Subject's Weigthed ImpGrid. It must be a "wimp" S3 object
#'        imported by  the \code{\link{importwimp}} function.
#' @param method Method for calculating centrality. You can use the simple
#'        method with "simple", normalized with "norm", weighted with "weigth",
#'        normalized weighted with "wnorm" and the ego density method with "ego".
#'        Default is Weigthed Method.
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
#' degree_index(example.wimp)
#'

degree_index <- function(wimp, method="weight"){

  lpoles <- wimp$constructs[[1]]
  rpoles <- wimp$constructs[[2]]
  poles <- wimp$constructs[[3]]

  wmat <- wimp$scores[[3]]
  N <- dim(wmat)[1]

  if(method == "simple" | method == "norm" | method == "ego"){
    wmat.1 <- wmat/wmat
    wmat.1[is.nan(wmat.1)] <- 0
    Cout <- rowSums(wmat.1)
    Cin <- colSums(wmat.1)
  }

  if(method == "weight" | method == "wnorm"){
    Cout <- rowSums(abs(wmat))
    Cin <- colSums(abs(wmat))
  }

  if(method == "norm" | method == "wnorm"){
    Cout <- Cout/(N-1)
    Cin <- Cin/(N-1)
  }

  if(method == "ego"){
    Cout <- Cout/(N*(N-1))
    Cin <- Cin/(N*(N-1))
  }

  names(Cout) <- poles
  names(Cin) <- poles

  result <- cbind(Cout, Cin , Cout + Cin)
  rownames(result) <- poles
  colnames(result) <- c("Out","In", "All")
  return(result)
}

# Distance Matrix ---------------------------------------------------------

#' Distance Matrix -- dismatrix()
#'
#' @description Function that calculates the shortest distance between each of
#'              the pairs of digraph constructions.
#'
#' @param wimp  Subject's Weigthed ImpGrid. It must be a "wimp" S3 object
#'        imported by the \code{\link{importwimp}} function.
#' @param mode Method to calculate the distances depending on the direction of
#'        the edges. With "out" we calculate them respecting the direction of the edges,
#'        "in" through the inverse of the direction of the edges and "all" without
#'        taking into account the direction. Default is "out"
#'
#' @author Alejandro Sanfeliciano
#'
#' @return Returns the digraph distance matrix. Matrix that contains the
#'         distances of the shortest paths from one construct to another.
#'
#' @export
#'

dismatrix <- function(wimp,mode="out"){

  poles <- wimp$constructs[[3]]
  wmat <- wimp$scores[[3]]

  G <- igraph::graph.adjacency(wmat,mode = "directed",weighted = T)

  result <- igraph::shortest.paths(G, weights = NA,mode = mode)

  rownames(result) <- poles
  colnames(result) <- poles

  return(result)
}

# Closeness Centrality Index ----------------------------------------------

#' Closeness index -- close_index()
#'
#' @description Function to calculate the closeness of a construct to the rest
#'              of the constructs within the digraph.
#'
#' @param wimp Subject's Weigthed ImpGrid. It must be a "wimp" S3 object
#'        imported by the \code{\link{importwimp}} function.
#' @param norm If TRUE, the values will be standardized. Default is TRUE.
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
#' close_index(example.wimp)
#' close_index(example.wimp, norm = FALSE)
#'

close_index <- function(wimp, norm = TRUE){

  lpoles <- wimp$constructs[[1]]
  rpoles <- wimp$constructs[[2]]
  poles <- wimp$constructs[[3]]

  dist <- dismatrix(wimp)
  N <- dim(dist)[1]

  result <- 1/(rowSums(dist))

  if(norm){
    result <- (N-1)/(rowSums(dist))
  }

  result <- matrix(result)
  rownames(result) <- poles
  colnames(result) <- "Closeness"

  return(result)
}

# Betweeness Centrality Index ---------------------------------------------

#' Betweeness index -- betw_index()
#'
#' @description Function that calculates the betweenness of each of the
#'              constructs. This is the number of times a geodesic path (shortest path)
#'              between two other constructs passes through that construct in the digraph.
#'
#' @param wimp  Subject's Weigthed ImpGrid. It must be a "wimp" S3 object
#'        imported by the \code{\link{importwimp}} function.
#' @param norm If TRUE, the values will be standardized. Default is TRUE.
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
#' betw_index(example.wimp)
#' betw_index(example.wimp, norm = FALSE)
#'

betw_index <- function(wimp,norm=TRUE){

  lpoles <- wimp$constructs[[1]]
  rpoles <- wimp$constructs[[2]]
  poles <- wimp$constructs[[3]]


  wmat <- wimp$scores[[3]]

  G <- igraph::graph.adjacency(wmat,mode = "directed",weighted = T)

  result <- igraph::betweenness(G,normalized = norm,weights = NA )

  result <- matrix(result)
  rownames(result) <- poles
  colnames(result) <- "Betweenness"

  return(result)
}

# PH Centrality Index ------------------------------------------------------------

#' Presence and Hierarchy Indices -- ph_index()
#'
#' @description This function computes the presence (P) and hierarchy (H) indices for constructs within an weigthed implication grid.
#'              These indices represent the frequency of occurrence and influence on other constructs, respectively.
#'              The function allows for different methods of standardization based on the context of the constructs.
#'
#' @param wimp A WIMP object containing an implication grid and associated constructs.
#' @param method A character string specifying the method used to calculate the degree indices.
#'        Default is "wnorm". Acceptable values include "wnorm", "simple", "weight", or any other method
#'        implemented in the 'degree_index' function.
#' @param std A character string indicating how to standardize the P and H indices. Available options are:
#'        - 'none': No standardization (default).
#'        - 'vertices': Standardizes by the maximum total degree, which is calculated based on the number of vertices.
#'        - 'edges': Standardizes by the total number of edges.
#'        - 'max_edges': Standardizes by the maximum number of outgoing edges from any single vertex.
#'        - 'density': Adjusts P and H by the density of the grid, which considers the total edges possible versus actual.
#'
#' @author Carlos Hurtado and Alejandro Sanfeliciano
#'
#' @return A matrix with two columns, 'p' for presence and 'h' for hierarchy, containing the indices for each construct.
#'         If standardization is applied, these values are modified according to the selected method.
#'
#' @export
#'
#' @examples
#'
#' ph_index(example.wimp)
#' ph_index(example.wimp, std = TRUE)
#' ph_index(example.wimp, method = "wnorm", std = FALSE)
#'


ph_index <- function(wimp, method = "wnorm", std = 'none'){

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

  # Standardization
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

# Eigen Indices ------------------------------------------------------------

#' Eigenvalue Centrality index -- eigen_index()
#'
#' @description This function calculates centrality scores for constructs within a `wimp` object,
#'              based on the eigenvalue decomposition of a specified adjacency matrix from the WIMP scores.
#'              It supports analyzing centrality using the 'direct', 'weights', or 'implications' matrices.
#'              The centrality calculation is performed over the specified number of eigenvectors.
#'
#' @param wimp wimp An object of class 'wimp' (weighted implications grid)
#' @param matrix A character string specifying which matrix to use for the centrality analysis. Accepted values are
#'        'direct', 'weights', or 'implications'. Default is 'implications'.
#' @param num.vectors An integer specifying the number of eigenvectors to use for computing centrality scores.
#'
#' @author Carlos Hurtado
#'
#' @return A dataframe containing the constructs' names and their respective centrality scores.
#'
#' @export
#'
#' @examples
#'
#' eigen_index(example.wimp)
#'

eigen_index <- function(wimp, matrix = "weights", num.vectors = 2) {
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
    stop("num.vectors exceeds the number of available eigenvectors.")
  }

  # Calculate centrality using the specified number of eigenvectors. This calculation takes the absolute value
  # of the real part of each eigenvector, squares it, and then multiplies it by the real part of the corresponding
  # eigenvalue to compute centrality scores.
  centralidad <- Reduce(`+`, lapply(1:num.vectors, function(i) {
    Re(results$vectors[, i])^2 * Re(results$values[i])
  }))

  # Create a dataframe with the centrality results
  df.centrality <- data.frame(
    Constructs = wimp$constructs$constructs,
    Eigenvalues = abs(centralidad)
  )

  return(df.centrality)
}

# PH Plot ------------------------------------------------

#' Scatter Plot of Constructs in PH Space -- ph_plot()
#'
#' @description This function generates a scatter plot of constructs in the P-H space,
#' where P represents Presence (frequency of the construct) and H represents Hierarchy
#' (influence of the construct).
#'
#' @param wimp An object of class 'wimp', which contains an implication grid
#'        and associated constructs.
#'
#' @param text.size Size of the text labels. Default is 1.
#'
#' @param ... Additional arguments are passed from \code{\link{ph_index}} function.
#'
#' @author Carlos Hurtado and Alejandro Sanfeliciano
#'
#' @return A Plotly object representing the generated scatter plot.
#'
#' @import plotly
#' @export
#'
#' @examples
#'
#' ph_plot(example.wimp)
#'

ph_plot <- function(wimp, text.size = 1, ...) {

  # Extract the Mahalobis distance matrix for the given wimp
  phm.mat <- ph_index(wimp)
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


  phm.mat.df$color <- .construct.colors(wimp = wimp, mode = "red/green")[,"color"]

  # Shapes for non-viable area if requested
  shapes <-
    list(
      list(type = "path", path = paste("M 0,0 L", limit, ",", limit, " L0,", limit, " Z"),
           fillcolor = "#FFD97D", opacity = 0.2, line = list(color = "#FA9D13")),
      list(type = "path", path = paste("M 0,0 L", limit, ",", -limit, " L0,", -limit, " Z"),
           fillcolor = "#FFD97D", opacity = 0.2, line = list(color = "#FA9D13")),
      list(type = "line", x0 = 0, y0 = 0, x1 = limit, y1 = limit,
           xref = "x", yref = "y", line = list(color = "#FFD97D", width = 1, dash = "dash")),
      list(type = "line", x0 = 0, y0 = 0, x1 = limit, y1 = -limit,
           xref = "x", yref = "y", line = list(color = "#FFD97D", width = 1, dash = "dash"))
    )

  # Initialize Plotly graph
  p <- plot_ly()

  # Set the layout of the graph
  p <- p %>%
    layout(title = '',
           xaxis = list(title = 'PRESENCE'),
           yaxis = list(title = 'HIERARCHY'),
           plot_bgcolor = "white",
           font = list(family = "Arial"),
           showlegend = FALSE,
           shapes = shapes)
  p <- p %>%
      add_annotations(data = phm.mat.df, x = ~p, y = ~h, text = ~self.constr,
                      hovertext = ~paste('Constructo:', constructo, '\nP:', p, 'H:', h), hoverinfo = 'text',
                      font = list(size = 12 * text.size, color = 'black'),
                      showarrow = FALSE, xanchor = 'center', yanchor = 'bottom',
                      yshift = 5)
  p <- p %>%
    add_markers(data = phm.mat.df, x = ~p, y = ~h,
                marker = list(color = phm.mat.df$color, size = 7, line = list(color = 'black', width = 1)),
                text = ~paste('P:', p, '; H:', h), hoverinfo = 'text')

  return(p)
}
