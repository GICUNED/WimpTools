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

ph_index <- function(wimp, method = "weight", std = FALSE){

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

  # Standardization based on number of graph edges
  if (std){
    #max.total.deg <- max(c.io[, 3])
    edges <- sum(c.io[, 2])
    ph.mat <- ph.mat / edges
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
#'   the PH matrix. Default is `"weight"`.
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

mahalanobis_index <- function(wimp, method = "weight", std = FALSE, sign.level = 0.2){

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

  # Annotate observations as "central" contructs based on the chi-square cutoff
  central <- ifelse(phm.mat[,"m.dist"] > chi.square.cutoff, TRUE, FALSE)

  # Add column to the results matrix
  phmc.mat <- cbind(phm.mat, central)
  return(phmc.mat)
}

# Graphically represent constructs ------------------------------------------------

#' Scatter Plot of Constructs in P-H Space
#'
#' This function generates a scatter plot of constructs in the P-H space,
#' where P represents Presence (frequency of the construct) and H represents Hierarchy
#' (influence of the construct). Central constructs are highlighted in red color, and
#' peripheral constructs in another, facilitating their visual identification.
#'
#' @param phm.mat A matrix where each row represents a construct and contains
#'        the P and H coordinates of the construct, as well as an indication of whether the
#'        construct is central (1) or not (0). The matrix must have row names,
#'        which are used to label the constructs in the graph.
#'
#' @param mark.nva Boolean value that specifies if non-viable areas are to be marked in the
#'        graphic. Non-viable areas are parts of the P-H space where constructs cannot logically exist.
#'
#' @param mark.cnt Boolean value that specifies if central constructs have to be highlighted.
#'        Central constructs are depicted with a distinct color and symbol to differentiate them from
#'        peripheral constructs.
#'
#' @param show.points Boolean value that specifies whether points should be displayed or not
#'
#' @return A Plotly object representing the generated scatter plot.
#'
#' @export
#'
#' @examples
#' graph_ph(phm.mat, mark.nva = TRUE, mark.cnt = TRUE)

graph_ph <- function(phm.mat, mark.nva = TRUE, mark.cnt = TRUE, show.points = TRUE) {

  # Convert the matrix to a dataframe
  phm.mat.df <- as.data.frame(phm.mat)
  # Assign the names of constructs from the row names of the matrix
  phm.mat.df$constructo <- rownames(phm.mat)
  # Colors for the constructs (default: peripheral constructs are non-central)
  colors <- viridis::viridis(n = nrow(phm.mat.df), option = "viridis")
  # Limits for the graph by the largest value of P or H dimensions. We add a small margin
  limit <- max(abs(phm.mat.df$p), abs(phm.mat.df$h)) * 1.1

  # Round P and H values
  phm.mat.df$p <- round(phm.mat.df$p, 3)
  phm.mat.df$h <- round(phm.mat.df$h, 3)

  # Add a new column for label color
  if (mark.cnt)
    phm.mat.df$label.color <- ifelse(phm.mat.df$central == 1, 'red', 'black')
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
    if (mark.cnt) {
      p <- p %>%
        add_markers(data = phm.mat.df[phm.mat.df$central == 0, ], x = ~p, y = ~h,
                    marker = list(color = colors, size = 10, line = list(color = 'black', width = 1)),
                    text = ~paste('P:', p, '; H:', h), hoverinfo = 'text') %>%
        add_markers(data = phm.mat.df[phm.mat.df$central == 1, ], x = ~p, y = ~h,
                    marker = list(color = 'orangered', size = 12, symbol = "circle-dot", line = list(color = 'black', width = 1)),
                    text = ~paste('P:', p, '; H:', h), hoverinfo = 'text')
    } else {
      p <- p %>%
        add_markers(data = phm.mat.df, x = ~p, y = ~h,
                    marker = list(color = colors, size = 10, line = list(color = 'black', width = 1)),
                    text = ~paste('P:', p, '; H:', h), hoverinfo = 'text')
    }
  }

  # Add annotations (labels) for each point
  #p <- p %>% add_annotations(data = phm.mat.df, x = ~p, y = ~h, text = ~constructo,
  #                           font = list(size = 12, color = ~label.color),
  #                           showarrow = FALSE, xanchor = 'center', yanchor = 'bottom')

  p <- p %>%
    add_annotations(
      data = phm.mat.df, x = ~p, y = ~h, text = ~constructo,
      hovertext = ~paste('Constructo:', constructo, '\nP:', p, 'H:', h), hoverinfo = 'text',
      #font = list(size = 12, color = ifelse(phm.mat.df$central == 1 & mark.cnt, 'red', 'black')),
      font = list(size = 12, color = 'black'),
      showarrow = FALSE, xanchor = 'center', yanchor = 'bottom'
    )

  return(p)
}
