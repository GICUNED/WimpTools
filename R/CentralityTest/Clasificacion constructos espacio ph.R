## CENTRALITY FUNCTIONS EDITADA PARA CLASIFICAR CONSTRUCTOS EN ESPACIO P_H ##

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

#MEDIANA +/- SD/2
ph_index_edicion2 <- function(wimp, method = "weight", std = 'none'){

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


  #data frame
  ph.mat.df <- as.data.frame(ph.mat)

  # Calcular la mediana y la desviación estándar de la columna H
  median_h <- median(ph.mat.df$h)
  std_dev_h <- sd(ph.mat.df$h)

  median_p <- median(ph.mat.df$p)
  std_dev_p <- sd(ph.mat.df$p)


  Descripcion <- ifelse(ph.mat.df$h > median(ph.mat.df$h) + (sd(ph.mat.df$h) / 2) &
            ph.mat.df$p > median(ph.mat.df$p) + (sd(ph.mat.df$p) / 2),
          "Nuclear, alta presencia",
            ifelse(ph.mat.df$h < median(ph.mat.df$h) - (sd(ph.mat.df$h) / 2) &
            ph.mat.df$p > median(ph.mat.df$p) + (sd(ph.mat.df$p) / 2),
          "Sintomático, alta presencia",
            ifelse(ph.mat.df$h <= median(ph.mat.df$h) + (sd(ph.mat.df$h) / 2) &
            ph.mat.df$h >= median(ph.mat.df$h) - (sd(ph.mat.df$h) / 2) &
            ph.mat.df$p > median(ph.mat.df$p) + (sd(ph.mat.df$p) / 2),
          "Mediador, alta presencia",
            ifelse(ph.mat.df$p <= median(ph.mat.df$p) + (sd(ph.mat.df$p) / 2) &
            ph.mat.df$p > median(ph.mat.df$p) - (sd(ph.mat.df$p) / 2) &
            ph.mat.df$h > median(ph.mat.df$h) + (sd(ph.mat.df$h) / 2),
          "Nuclear, moderada presencia",
            ifelse(ph.mat.df$p <= median(ph.mat.df$p) + (sd(ph.mat.df$p) / 2) &
            ph.mat.df$p > median(ph.mat.df$p) - (sd(ph.mat.df$p) / 2) &
            ph.mat.df$h < median(ph.mat.df$h) - (sd(ph.mat.df$h) / 2),
          "Sintomático, moderada presencia",
            ifelse(ph.mat.df$p <= median(ph.mat.df$p) + (sd(ph.mat.df$p) / 2) &
            ph.mat.df$p > median(ph.mat.df$p) - (sd(ph.mat.df$p) / 2) &
            ph.mat.df$h >= median(ph.mat.df$h) - (sd(ph.mat.df$h) / 2) &
            ph.mat.df$h <= median(ph.mat.df$h) + (sd(ph.mat.df$h) / 2),
          "Mediador, moderada presencia",
          "Superficial"
          ))))))

    ph.mat.df2 <- cbind(ph.mat.df, Descripcion)

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

  return(ph.mat.df2)
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

mahalanobis_index_edicion2  <- function(wimp, method = "weight", std = 'none', sign.level = 0.2){

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

  # Annotate observations as "central" constructs based on the chi-square cutoff
  # and being "to the right" on the P axis
  central <- ifelse(phm.mat[,"m.dist"] > chi.square.cutoff & phm.mat[,"p"] > mean(phm.mat[,"p"]), TRUE, FALSE)

  # Add column to the results matrix
  phmc.mat <- cbind(phm.mat, central)
  return(phmc.mat)
}

# Graphically represent constructs ------------------------------------------------

#' Scatter Plot of Wimps Constructs in P-H Space
#'
#' This function generates a scatter plot of constructs of a given wimp in the P-H space,
#' where P represents Presence (frequency of the construct) and H represents Hierarchy
#' (influence of the construct). Central constructs are highlighted in red color, and
#' peripheral constructs in another, facilitating their visual identification.
#'
#' @param wimp A `wimp` object containing the data from which the PH matrix
#'   is derived.
#'
#' @param method A character string specifying the method used for calculating
#'   the PH matrix. Default is `"weight"`.
#'
#' @param std A logical value indicating whether the data should be
#'   standardized before calculating the Mahalanobis distance. Default is `FALSE`.
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
#'
graph_ph_edicion2  <- function(..., mark.nva = TRUE, mark.cnt = TRUE, show.points = TRUE) {

  # Extract the Mahalobis distance matrix for the given wimp
  phm.mat <- mahalanobis_index(...)
  # Convert the matrix to a dataframe
  phm.mat.df <- as.data.frame(phm.mat)
  # Assign the names of constructs from the row names of the matrix
  phm.mat.df$constructo <- rownames(phm.mat)
  # Limits for the graph by the largest value of P or H dimensions. We add a small margin
  limit <- max(abs(phm.mat.df$p), abs(phm.mat.df$h)) * 1.1

  # Round P and H values
  phm.mat.df$p <- round(phm.mat.df$p, 3)
  phm.mat.df$h <- round(phm.mat.df$h, 3)


  # Calcular la mediana y la desviación estándar de la columna H
  median_h <- median(phm.mat.df$h)
  std_dev_h <- sd(phm.mat.df$h)

  median_p <- median(phm.mat.df$p)
  std_dev_p <- sd(phm.mat.df$p)


  # Calcular el límite superior usando la mediana y la desviación estándar / 2 (la mitad de las desviación)
  upper_limit <- median_h+(std_dev_h/2)
  lower_limit <- median_h-(std_dev_h/2)

  right_limit <- median_p+(std_dev_p/2)
  left_limit <- median_p-(std_dev_p/2)

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
           xref = "x", yref = "y", line = list(color = "darkgreen", width = 1, dash = "dash")),



        # Agregar una línea horizontal adicional upper_limit
    list(type = "line", x0 = 0, y0 = upper_limit, x1 = limit, y1 = upper_limit,
         xref = "x", yref = "y", line = list(color = "blue", width = 1, dash = "dot")),

    # Agregar una línea horizontal adicional lower_limit

    list(type = "line", x0 = 0, y0 = lower_limit, x1 = limit, y1 = lower_limit,
         xref = "x", yref = "y", line = list(color = "blue", width = 1, dash = "dot")),

    # Agregar una línea vertical que divide el espacio en dos partes a partir de la mediana de P right_limit
    list(type = "line", x0 = right_limit, y0 = -limit, x1 = right_limit, y1 = limit,
         xref = "x", yref = "y", line = list(color = "blue", width = 1, dash = "dot")),

    # Agregar una línea vertical que divide el espacio en dos partes a partir de la mediana de P left_limit
    list(type = "line", x0 = left_limit, y0 = -limit, x1 = left_limit, y1 = limit,
         xref = "x", yref = "y", line = list(color = "blue", width = 1, dash = "dot")
    ))
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

    if (mark.cnt) {
      p <- p %>%
        add_markers(data = phm.mat.df, x = ~p, y = ~h,
                    marker = list(color = ~color, size = 10,
                                  symbol = ifelse(phm.mat.df$central == 1, "star", "circle"),
                                  line = list(color = 'black', width = ifelse(phm.mat.df$central == 1, 2, 1))),
                    text = ~paste('P:', p, '; H:', h), hoverinfo = 'text')
    } else {
      p <- p %>%
        add_markers(data = phm.mat.df, x = ~p, y = ~h,
                    #marker = list(color = colors, size = 10, line = list(color = 'black', width = 1)),
                    marker = list(color = phm.mat.df$color, size = 10, line = list(color = 'black', width = 1)),
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
      font = list(size = 12, color = phm.mat.df$label.color),
      #font = list(size = 12, color = .label.color(phm.mat.df$central == 1 & mark.cnt == TRUE)),
      showarrow = FALSE, xanchor = 'center', yanchor = 'bottom'
    )

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

#MEDIA +/- SD/2
ph_index_edicion3 <- function(wimp, method = "weight", std = 'none'){

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


  #data frame
  ph.mat.df <- as.data.frame(ph.mat)

  # Calcular la mediana y la desviación estándar de la columna H
  mean_h <- mean(ph.mat.df$h)
  std_dev_h <- sd(ph.mat.df$h)

  mean_p <- mean(ph.mat.df$p)
  std_dev_p <- sd(ph.mat.df$p)


  Descripcion <- ifelse(ph.mat.df$h > mean(ph.mat.df$h) + (sd(ph.mat.df$h) / 2) &
                          ph.mat.df$p > mean(ph.mat.df$p) + (sd(ph.mat.df$p) / 2),
                      "Nuclear, alta presencia",
                        ifelse(ph.mat.df$h < mean(ph.mat.df$h) - (sd(ph.mat.df$h) / 2) &
                        ph.mat.df$p > mean(ph.mat.df$p) + (sd(ph.mat.df$p) / 2),
                      "Sintomático, alta presencia",
                        ifelse(ph.mat.df$h <= mean(ph.mat.df$h) + (sd(ph.mat.df$h) / 2) &
                        ph.mat.df$h >= mean(ph.mat.df$h) - (sd(ph.mat.df$h) / 2) &
                        ph.mat.df$p > mean(ph.mat.df$p) + (sd(ph.mat.df$p) / 2),
                      "Mediador, alta presencia",
                       ifelse(ph.mat.df$p <= mean(ph.mat.df$p) + (sd(ph.mat.df$p) / 2) &
                       ph.mat.df$p > mean(ph.mat.df$p) - (sd(ph.mat.df$p) / 2) &
                       ph.mat.df$h > mean(ph.mat.df$h) + (sd(ph.mat.df$h) / 2),
                      "Nuclear, moderada presencia",
                       ifelse(ph.mat.df$p <= mean(ph.mat.df$p) + (sd(ph.mat.df$p) / 2) &
                       ph.mat.df$p > mean(ph.mat.df$p) - (sd(ph.mat.df$p) / 2) &
                       ph.mat.df$h < mean(ph.mat.df$h) - (sd(ph.mat.df$h) / 2),
                     "Sintomático, moderada presencia",
                       ifelse(ph.mat.df$p <= mean(ph.mat.df$p) + (sd(ph.mat.df$p) / 2) &
                       ph.mat.df$p > mean(ph.mat.df$p) - (sd(ph.mat.df$p) / 2) &
                       ph.mat.df$h >= mean(ph.mat.df$h) - (sd(ph.mat.df$h) / 2) &
                       ph.mat.df$h <= mean(ph.mat.df$h) + (sd(ph.mat.df$h) / 2),
                      "Mediador, moderada presencia",
                      "Superficial"
                      ))))))

  ph.mat.df2 <- cbind(ph.mat.df, Descripcion)

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

  return(ph.mat.df2)
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

mahalanobis_index_edicion3  <- function(wimp, method = "weight", std = 'none', sign.level = 0.2){

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

  # Annotate observations as "central" constructs based on the chi-square cutoff
  # and being "to the right" on the P axis
  central <- ifelse(phm.mat[,"m.dist"] > chi.square.cutoff & phm.mat[,"p"] > mean(phm.mat[,"p"]), TRUE, FALSE)

  # Add column to the results matrix
  phmc.mat <- cbind(phm.mat, central)
  return(phmc.mat)
}

# Graphically represent constructs ------------------------------------------------

#' Scatter Plot of Wimps Constructs in P-H Space
#'
#' This function generates a scatter plot of constructs of a given wimp in the P-H space,
#' where P represents Presence (frequency of the construct) and H represents Hierarchy
#' (influence of the construct). Central constructs are highlighted in red color, and
#' peripheral constructs in another, facilitating their visual identification.
#'
#' @param wimp A `wimp` object containing the data from which the PH matrix
#'   is derived.
#'
#' @param method A character string specifying the method used for calculating
#'   the PH matrix. Default is `"weight"`.
#'
#' @param std A logical value indicating whether the data should be
#'   standardized before calculating the Mahalanobis distance. Default is `FALSE`.
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
#'
graph_ph_edicion3  <- function(..., mark.nva = TRUE, mark.cnt = TRUE, show.points = TRUE) {

  # Extract the Mahalobis distance matrix for the given wimp
  phm.mat <- mahalanobis_index(...)
  # Convert the matrix to a dataframe
  phm.mat.df <- as.data.frame(phm.mat)
  # Assign the names of constructs from the row names of the matrix
  phm.mat.df$constructo <- rownames(phm.mat)
  # Limits for the graph by the largest value of P or H dimensions. We add a small margin
  limit <- max(abs(phm.mat.df$p), abs(phm.mat.df$h)) * 1.1

  # Round P and H values
  phm.mat.df$p <- round(phm.mat.df$p, 3)
  phm.mat.df$h <- round(phm.mat.df$h, 3)


  # Calcular la mediana y la desviación estándar de la columna H
  mean_h <- mean(phm.mat.df$h)
  std_dev_h <- sd(phm.mat.df$h)

  mean_p <- mean(phm.mat.df$p)
  std_dev_p <- sd(phm.mat.df$p)


  # Calcular el límite superior usando la mediana y la desviación estándar / 2 (la mitad de las desviación)
  upper_limit <- mean_h+(std_dev_h/2)
  lower_limit <- mean_p-(std_dev_h/2)

  right_limit <- mean_h+(std_dev_p/2)
  left_limit <- mean_p-(std_dev_p/2)

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
           xref = "x", yref = "y", line = list(color = "darkgreen", width = 1, dash = "dash")),



      # Agregar una línea horizontal adicional upper_limit media de h
      list(type = "line", x0 = 0, y0 = upper_limit, x1 = limit, y1 = upper_limit,
           xref = "x", yref = "y", line = list(color = "blue", width = 1, dash = "dot")),

      # Agregar una línea horizontal adicional lower_limit media de h

      list(type = "line", x0 = 0, y0 = lower_limit, x1 = limit, y1 = lower_limit,
           xref = "x", yref = "y", line = list(color = "blue", width = 1, dash = "dot")),

      # Agregar una línea vertical que divide el espacio en dos partes a partir de la media de P right_limit
      list(type = "line", x0 = right_limit, y0 = -limit, x1 = right_limit, y1 = limit,
           xref = "x", yref = "y", line = list(color = "blue", width = 1, dash = "dot")),

      # Agregar una línea vertical que divide el espacio en dos partes a partir de la media de P left_limit
      list(type = "line", x0 = left_limit, y0 = -limit, x1 = left_limit, y1 = limit,
           xref = "x", yref = "y", line = list(color = "blue", width = 1, dash = "dot")
      ))
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

    if (mark.cnt) {
      p <- p %>%
        add_markers(data = phm.mat.df, x = ~p, y = ~h,
                    marker = list(color = ~color, size = 10,
                                  symbol = ifelse(phm.mat.df$central == 1, "star", "circle"),
                                  line = list(color = 'black', width = ifelse(phm.mat.df$central == 1, 2, 1))),
                    text = ~paste('P:', p, '; H:', h), hoverinfo = 'text')
    } else {
      p <- p %>%
        add_markers(data = phm.mat.df, x = ~p, y = ~h,
                    #marker = list(color = colors, size = 10, line = list(color = 'black', width = 1)),
                    marker = list(color = phm.mat.df$color, size = 10, line = list(color = 'black', width = 1)),
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
      font = list(size = 12, color = phm.mat.df$label.color),
      #font = list(size = 12, color = .label.color(phm.mat.df$central == 1 & mark.cnt == TRUE)),
      showarrow = FALSE, xanchor = 'center', yanchor = 'bottom'
    )

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

construct_colors_edicion3 <- function(wimp, mode){
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



