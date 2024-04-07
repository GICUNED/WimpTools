
## EMOTIONAL ADJUSTMENT ##

# calculate_adjustment_init_end ------------------------------------------------------------
#'
#' This function calculates the Euclidean distance and the correlation between the end vector and init vector.
#'
#' @param vector_init Object containing the initial vector.
#' @param vector_end Object containing the end vector.
#' @param normalize Logical. If TRUE, normalize the distance; otherwise, use the raw Euclidean distance.
#' @param MaxScale Number that contains the maximum value of each dimension of the vector. Used to calculate the normalized Euclidean distance between vector_init and vector_end
#' @param filtered_constructs In case we want to use only specific values (e.g use nuclear constructs only). Vector with 0s and 1s.
#'Only the values from the 'vector_end' and 'vector_init' vectors at positions where there is a 1 will be included.
#'By default, the 'vector_end' and 'vector_init' vectors are kept.#'
#' @return A list with the calculated values: distance, correlation and magnitude (as the adjustment index).
#' @examples calculate_adjustment_init_end (vector_init, vector_end,,, c(1,1,0,0,1,0,1,0))
#' @export

calculate_adjustment_init_end <- function(vector_init, vector_end, normalize = TRUE, MaxScale = 1, filtered_constructs = c(1)) {

  # Use only filtered_constructs:
  if (length(filtered_constructs) != 0) {
    n <- length(filtered_constructs)
    elems <- 1:n

    for (i in elems) {
      if (filtered_constructs[n+1-i] == 0) {
        vector_end <- vector_end[-n-1+i]
        vector_init <- vector_init[-n-1+i]
      }
    }
  }

  # Calculate the euclidean distance between both vectors
  eu_distance <- sqrt(sum((vector_end - vector_init)^2))

  # Normalize the Euclidean distance if specified
  distance <- if (normalize) eu_distance / (2 * MaxScale * sqrt(length(vector_end)))  else eu_distance

  # Calculate the correlation between the vectors (cosine)
  #WARNING: Return NaN if one of the vectors is "0". This may happen with the filtering.
  correlation <- sum(vector_end * vector_init) / (norm(vector_end, type = "2") * norm(vector_init, type = "2"))

  # General adjustment value (magnitude of the calculated vector)
  magnitude <- if (normalize) 1 - (sqrt((1-correlation)^2 + distance^2)) / sqrt(5) else NULL

  # Return the results
  return(list(distance = distance, correlation = correlation, magnitude=magnitude))
}

# calculate_adjustment_vector ------------------------------------------------------------
#'
#' This function calculates the normalized length of a vector.
#'
#' @param vector Object containing the vector.
#' @param normalize Logical. If TRUE, normalize the length; otherwise, use the raw length.
#' @param MaxScale Number that contains the maximum value of each dimension of the vector. Used to calculate the normalized length of the vector.
#' @param filtered_constructs In case we want to use only specific values (e.g use nuclear constructs only). Vector with 0s and 1s.
#'Only the values from the vector at positions where there is a 1 will be included.
#' @return The length (normalized or not) of the vector.
#' @examples calculate_adjustment_vector (vector,,, c(1,1,0,0,1,0,1,0))
#' @export

calculate_adjustment_vector <- function(vector, normalize = TRUE, MaxScale = 1, filtered_constructs = c(1)) {

  # Use only filtered_constructs:
  if (length(filtered_constructs) != 0) {
    n <- length(filtered_constructs)
    elems <- 1:n

    for (i in elems) {
      if (filtered_constructs[n+1-i] == 0) {
        vector <- vector[-n-1+i]
      }
    }
  }

  # Calculate the length of the vector (Euclidean norm)
  eu_length <- norm(vector, type = "2")

  # Normalize the Euclidean distance if specified
  eu_length <- if (normalize) eu_length / (MaxScale * sqrt(length(vector)))  else eu_length

  # Return the result
  return(eu_length)
}



## EMOTIONAL ADJUSTMENT ##

# calculate_adjustment_wimp ------------------------------------------------------------
#'
#' This function calculates the Euclidean distance, the correlation and the adjustment index for:
#' the self,
#' the hypothetical self when each construct move from self to ideal (considering and not considering the original change in the construct),
#' the hypothetical self when each construct move one step from self to ideal (considering and not considering the original change in the construct),
#' all above against the ideal.
#' Also the adjustment index as the indefinition of self and ideal.
#' @param wimp Object containing ideal, self vectors and weights matrix.
#' @param filtered_constructs In case we want to use only specific values (e.g use nuclear constructs only). Vector with 0s and 1s.
#'Only the values from the vectors at positions where there is a 1 will be included.
#' @return A list with the calculated values: distance, correlation and magnitude, for self and every hypothetical self, and the adjustment index as the indefinition of self and ideal.
#' @examples calculate_adjustment_wimp (wimp,c(1,1,0,0,1,0,1,0))
#' @export
calculate_adjustment_wimp <- function(wimp, filtered_constructs = c(1)) {
  adjustment <- list()

  # Get the ideal and self vectors
  vector_ideal <- wimp$ideal[[2]]
  vector_self <- wimp$self[[2]]

  # CALCULATE ADJUSTMENT OF SELF
  adjustment$self <- list()
  adjustment$self <- calculate_adjustment_init_end(vector_self, vector_ideal, TRUE, 1, filtered_constructs = filtered_constructs)

  # Add labels to adjustment$self
  names(adjustment$self) <- c("distance", "correlation", "adjustment")

  # CALCULATE ADJUSTMENT OF HYPOTHETICAL SELF. 4 options:
  # ToIdeal: Contains the adjustment index for all the hypothetical-self, when each construct change to the ideal
  # ToIdeal2: Same as ToIdeal, but removing the change in the construct
  # OneStep: Contains the adjustment index for all the hypothetical-self, when each construct change 1-step to the ideal
  # OneStep2: Same as OneStep, but removing the change in the construct
  adjustment$future <- list(
    ToIdeal = list(),
    ToIdeal2 = list(),
    OneStep = list(),
    OneStep2 = list()
  )

  # Calculate vector_future
  # Calculate the difference between the ideal and self vectors
  difference <- vector_ideal - vector_self
  difference_OneStep <- (vector_ideal - vector_self)/abs(vector_ideal - vector_self)

  # Get the length of the difference vector
  length <- length(difference)

  # Build the matrix with the difference on the diagonal and zeros elsewhere
  change_vectors_ToIdeal <- diag(difference, nrow = length, ncol = length)
  # change_vectors_OneStep <- change_vectors_ToIdeal / abs(change_vectors_ToIdeal)
  change_vectors_OneStep <- diag(difference_OneStep, nrow = length, ncol = length) * 2 / 3

  # Get the weight matrix
  weight_matrix <- wimp$scores$weights

  # Multiply change_vectors by the weight matrix for the 2 options (ToIdeal and OneStep)
  result_matrix_ToIdeal <- change_vectors_ToIdeal %*% weight_matrix
  result_matrix_OneStep <- change_vectors_OneStep %*% weight_matrix

  # Compose vector_self to sum by rows
  vector_self_row_sum <- matrix(vector_self, nrow = length, ncol = length, byrow = TRUE)

  # Sum the result with the vector_self and change_vectors for all the options
  vector_future_ToIdeal <- vector_self_row_sum + change_vectors_ToIdeal + result_matrix_ToIdeal
  vector_future_ToIdeal2 <- vector_self_row_sum + result_matrix_ToIdeal
  vector_future_OneStep <- vector_self_row_sum + change_vectors_OneStep + result_matrix_OneStep
  vector_future_OneStep2 <- vector_self_row_sum + result_matrix_OneStep

  # Limit values in vector_future to be within -1 and +1 for all the options
  vector_future_ToIdeal <- pmax(pmin(vector_future_ToIdeal, 1), -1)
  vector_future_ToIdeal2 <- pmax(pmin(vector_future_ToIdeal2, 1), -1)
  vector_future_OneStep <- pmax(pmin(vector_future_OneStep, 1), -1)
  vector_future_OneStep2 <- pmax(pmin(vector_future_OneStep2, 1), -1)

  # Iterate over vectors in vector_future for all the options
  for (i in 1:nrow( vector_future_ToIdeal)) {
    # Calculate adjustment for vector_future and ideal vectors
    adjustment$future$ToIdeal[[i]] <- calculate_adjustment_init_end(vector_future_ToIdeal[i,], vector_ideal, TRUE, 1, filtered_constructs = filtered_constructs)
    # Add labels to adjustment$future
    names(adjustment$future$ToIdeal[[i]]) <- c("distance", "correlation", "adjustment")
  }
  for (i in 1:nrow( vector_future_ToIdeal2)) {
    # Calculate adjustment for vector_future and ideal vectors
    adjustment$future$ToIdeal2[[i]] <- calculate_adjustment_init_end(vector_future_ToIdeal2[i,], vector_ideal, TRUE, 1, filtered_constructs = filtered_constructs)
    # Add labels to adjustment$future
    names(adjustment$future$ToIdeal2[[i]]) <- c("distance", "correlation", "adjustment")
  }
  for (i in 1:nrow( vector_future_OneStep)) {
    # Calculate adjustment for vector_future and ideal vectors
    adjustment$future$OneStep[[i]] <- calculate_adjustment_init_end(vector_future_OneStep[i,], vector_ideal, TRUE, 1, filtered_constructs = filtered_constructs)
    # Add labels to adjustment$future
    names(adjustment$future$OneStep[[i]]) <- c("distance", "correlation", "adjustment")
  }
  for (i in 1:nrow( vector_future_OneStep2)) {
    # Calculate adjustment for vector_future and ideal vectors
    adjustment$future$OneStep2[[i]] <- calculate_adjustment_init_end(vector_future_OneStep2[i,], vector_ideal, TRUE, 1, filtered_constructs = filtered_constructs)
    # Add labels to adjustment$future
    names(adjustment$future$OneStep2[[i]]) <- c("distance", "correlation", "adjustment")
  }

  # CALCULATE ADJUSTMENT AS INDEFINITION OF SELF AND IDEAL

  adjustment$indefinition <- list(
    Self = calculate_adjustment_vector(vector_self, TRUE, 1, filtered_constructs = filtered_constructs),
    Ideal = calculate_adjustment_vector(vector_ideal, TRUE, 1, filtered_constructs = filtered_constructs)
  )
  # Add labels to adjustment$indefinition
  names(adjustment$indefinition) <- c("Indefinition of Self", "Indefinition of Ideal")

  # Return adjustment
  return(adjustment)

}



## EMOTIONAL ADJUSTMENT ##

# calculate_adjustment_self_ideal ------------------------------------------------------------
#'
#' This function calculates the Ecuclidean distance and the correlation between the ideal vector and self vector.
#'
#' @param wimp Object containing ideal and self vectors.
#' @param normalize Logical. If TRUE, normalize the distance; otherwise, use the raw Euclidean distance.
#' @param filtered_constructs In case we want to use only specific values (e.g use nuclear constructs only). Vector with 0s and 1s.
#'Only the values from the 'vector_ideal' and 'vector_self' vectors at positions where there is a 1 will be included.
#'By default, the 'vector_ideal' and 'vector_self' vectors are kept.#'
#' @return A list with the calculated values: distance and correlation.
#' @examples calculate_adjustment_self_ideal (wimp,TRUE,c(1,1,0,0,1,0,1,0))
#' @export

calculate_adjustment_self_ideal <- function(wimp, normalize = TRUE, filtered_constructs = c(1)) {

  # Get ideal and self vectors from the wimp variable
  vector_ideal <- wimp$ideal[[2]]
  vector_self <- wimp$self[[2]]


  # Use only filtered_constructs:
  if (length(filtered_constructs) != 0) {
    n <- length(filtered_constructs)
    elems <- 1:n

    for (i in elems) {
      if (filtered_constructs[n+1-i] == 0) {
        vector_ideal <- vector_ideal[-n-1+i]
        vector_self <- vector_self[-n-1+i]
      }
    }
  }

  # Calculate the euclidean distance between both vectors
  eu_distance <- sqrt(sum((vector_ideal - vector_self)^2))

  # Normalize the Euclidean distance if specified
  distance <- if (normalize) eu_distance / (2 * sqrt(length(vector_ideal)))  else eu_distance

  # Calculate the correlation between the vectors (cosine)
  correlation <- sum(vector_ideal * vector_self) / (norm(vector_ideal, type = "2") * norm(vector_self, type = "2"))

  # General adjustment value (magnitude of the calculated vector)
  magnitude <- if (normalize) 1 - (sqrt((1-correlation)^2 + distance^2)) / sqrt(5) else NULL

  # Return the results
  return(list(distance = distance, correlation = correlation, magnitude=magnitude))
}

#' plot_adjustment_self_ideal -----------------------------------------------------------------
#'
#' Plot graphs for ideal, self, and adjustment vectors.
#'
#' @param wimp Object containing ideal and self vectors.
#' @param calculated_values Results from the calculate_adjustment_self_ideal function. If empty get the wimp ones.
#' @param filtered_constructs In case we want to use only specific values (e.g use nuclear constructs only). The function is
#'        set to get a list with two vectors with the standardized ideal and self values "vector_ideal" and "vector_self".
#' @return The graphic vector for further use.
#'

plot_adjustment_self_ideal <- function(wimp, calculated_values = list(), filtered_constructs = c(1)) {

  # Get ideal and self vectors

  vector_ideal <- wimp$ideal[[2]]
  vector_self <- wimp$self[[2]]


  # Use only filtered_constructs:
  if (length(filtered_constructs) != 0) {
    n <- length(filtered_constructs)
    elems <- 1:n

    for (i in elems) {
      if (filtered_constructs[n + 1 - i] == 0) {
        vector_ideal <- vector_ideal[-n - 1 + i]
        vector_self <- vector_self[-n - 1 + i]
      }
    }
  }

  # Get the calculated graphic vector
  calculated_values <- if (length(calculated_values) == 0) calculate_adjustment_self_ideal(wimp) else calculated_values
  graphic_vector <- c(1, 0, calculated_values$correlation, calculated_values$distance)

  # Set ideal and self vectors values for the plot
  vector_ideal_norm <- norm(vector_ideal, type="2")
  vector_self_norm <- norm(vector_self, type="2")
  n <- length (vector_ideal)
  vector_ideal_x <- 0
  vector_ideal_y <- vector_ideal_norm / sqrt(n)
  theta <- acos(calculated_values$correlation) #radians
  vector_self_x <- vector_self_norm * sin(theta) / sqrt(n)
  vector_self_y <- vector_self_norm * calculated_values$correlation / sqrt(n)

  #
  # Set graphic for ideal and self vector
  #

  result <- c(45, 45, 45, 45, 45, 45, 45, 45)
  colors <- c("palegreen", "yellow","orange", "red", "red", "orange","yellow", "palegreen" )

  # if labels and no legend
  #alabels <- c("muy buen ajuste","buen ajuste","mal ajuste","muy mal ajuste",
  #            "muy mal ajuste","mal ajuste","buen ajuste","muy buen ajuste")
  #pie(result, main="Nivel de ajuste I", init.angle = 90, radius = 1, col=colors, labels=alabels)

  # if legend and no labels:
  pie(result, main="Nivel de ajuste", init.angle = 90, radius = 1, col=colors, labels=c(""))
  # draw the legend
  #legend(-2, 1, c("muy buen ajuste", "buen ajuste", "mal ajuste","muy mal ajuste"), fill=c("green", "yellow","orange", "red"))

  #circumferences aroud ideal vector
  xcenter <- 0
  ycenter <- vector_ideal_y
  theta <- seq(0, 2 * pi, length = 200)

  rad     <- vector_ideal_y + 1
  polygon(x=rad * cos(theta) + xcenter,
          y=rad * sin(theta) + ycenter,
          lwd=3, lty="solid", border='firebrick')

  rad     <- 0.75 * (vector_ideal_y + 1)
  polygon(x=rad * cos(theta) + xcenter,
          y=rad * sin(theta) + ycenter,
          lwd=3, lty="dashed", border='darkorange1')

  rad     <- 0.5 * (vector_ideal_y + 1)
  polygon(x=rad * cos(theta) + xcenter,
          y=rad * sin(theta) + ycenter,
          lwd=3, lty="dotted", border='khaki4')

  rad     <- 0.25 * (vector_ideal_y + 1)
  polygon(x=rad * cos(theta) + xcenter,
          y=rad * sin(theta) + ycenter,
          lwd=3, lty="twodash", border='darkgreen')

  # Draw ideal and self vectors
  arrows(0, 0, vector_ideal_x, vector_ideal_y, col = "blue", length = 0.1, lwd = 4)
  arrows(0, 0, vector_self_x, vector_self_y, col = "purple", length = 0.1,lwd = 4)
  text(-0, 0.5, expression(bold("Vector Ideal")), col = "blue")
  text(vector_self_x + 0.2, vector_self_y + 0.2, expression(bold("Vector Self")), col = "purple")

  #
  #
  # Set graphic for adjustment vector: origin (1, 0) and end point (correlation, distance)
  #
  plot(c(0, 1.5), c(0, 1.5), type = "n", xlab = "correlación", ylab = "distancia normalizada", asp = 1, xlim = c(-1, 1), ylim = c(0, 1.5))

  # Draw the different adjustement areas
  rect(-1, 0, 1, 1, density=20, col='red')
  rect(-0.5, 0, 1, 1, density=20, col='orange')
  rect(0.25, 0, 1, 0.75, density=20, col='yellow')
  rect(0.75, 0, 1, 0.5, density=20, col='green')

  # Define the margins of the graphic and the Abscissa and ordinate axis
  abline(v=c(-1,0,1), lwd = 2, lty=c(2,1,2), col=c("red","black","red"))
  abline(h=c(0,1), lwd = 2, lty=c(1,2), col=c("black","red"))

  # Draw the final vector [(1,0), (correlation, stand. distance)]
  tolerance <- 1e-10  # Define a small tolerance
  if (abs(calculated_values$correlation - 1) < tolerance && calculated_values$distance==0){
    points(1, 0, col = "blue", pch = 16, cex = 1.2)}
  else{
    arrows(1, 0, calculated_values$correlation, calculated_values$distance, col = "blue", length = 0.1)
  }

  text(1.2, 0.1, "Vector (1,0) a (correlación, distancia_normalizada)", col = "blue")

  # Return the graphic vector for further use if needed
  return(graphic_vector)

}

#' Adjustment_self_ideal ---------------------------------------------------------------------------------------
#'
#' This function calculates the Ecuclidean distance and the correlation between the ideal vector and self vector and plot corresponding graphs.
#'
#' @param wimp Object containing ideal and self vectors.
#' @param normalize Logical. If TRUE, normalize the distance; otherwise, use the raw Euclidean distance.
#' @return A list with the calculated values and the graphic vector.
#'

Adjustment_self_ideal <- function(wimp, normalize = TRUE, filtered_constructs = c(1)) {
  # Calculate adjustment and self-ideal values
  calculated_values <- calculate_adjustment_self_ideal(wimp, normalize = normalize, filtered_constructs = filtered_constructs )

  # Plot the graphs
  graphic_vector <- plot_adjustment_self_ideal(wimp, calculated_values, filtered_constructs = filtered_constructs)

  # Return data
  return(list(distance = calculated_values$distance, correlation = calculated_values$correlation, magnitude = calculated_values$magnitude, graphic = graphic_vector))
}


#' load_data_adjustment ------------------------------------------------------------
#'
#' This function extracts from 'data' the necessary values to study emotional adjustment.
#'
#' @return A data frame containing the extracted data for analysis. The data frame includes the following columns:
#'   ID - Identification code for each participant.
#'   dist.w1 - Distance measure for emotional adjustment in the WIMP test.
#'   corr.w1 - Correlation measure for emotional adjustment in the WIMP test.
#'   magn.w1 - Magnitude measure for emotional adjustment in the WIMP test.
#'   dist.w2 - Distance measure for emotional adjustment in the WIMP retest.
#'   corr.w2 - Correlation measure for emotional adjustment in the WIMP retest.
#'   magn.w2 - Magnitude measure for emotional adjustment in the WIMP retest.
#'   swls.testS - WLS test scores (Life Satisfaction).
#'   swls.retest - SWLS retest scores (Life Satisfaction).
#'   rse.testRSE -  test scores (Self-Esteem).
#'   rse.retestRSE - retest scores (Self-Esteem).
#'   tas.testTAS - test scores (Alexithymia).
#'   tas.retestTAS - retest scores (Alexithymia).
#'   mcq.testMCQ - test scores (Metacognition).
#'   mcq.retestMCQ - retest scores (Metacognition).

load_data_adjustment <- function () {
  data(data)
  df_data <- data.frame(
    ID = character(),
    dist.w1 = numeric(),
    corr.w1 = numeric(),
    magn.w1 = numeric(),
    dist.w2 = numeric(),
    corr.w2 = numeric(),
    magn.w2 = numeric(),
    swls.test = numeric(),
    swls.retest = numeric(),
    rse.test = numeric(),
    rse.retest = numeric(),
    tas.test = numeric(),
    tas.retest = numeric(),
    mcq.test = numeric(),
    mcq.retest = numeric()
  )

  for (i in 1:nrow(data$dataset)) {
    id_index <- match(data$dataset$ID[i], names(data$grids))
    if (!is.na(id_index)) {
      adjustw1 <- calculate_adjustment_self_ideal(data$grids[[id_index]][[3]])
      adjustw2 <- calculate_adjustment_self_ideal(data$grids[[id_index]][[4]])

      row_data <- data.frame(
        ID = data$dataset$ID[i],
        dist.w1 = adjustw1[1],
        corr.w1 = adjustw1[2],
        magn.w1 = adjustw1[3],
        dist.w2 = adjustw2[1],
        corr.w2 = adjustw2[2],
        magn.w2 = adjustw2[3],
        swls.test = data$dataset$swls.test[i],
        swls.retest = data$dataset$swls.retest[i],
        rse.test = data$dataset$rse.test[i],
        rse.retest = data$dataset$rse.retest[i],
        tas.test = data$dataset$tas.test[i],
        tas.retest = data$dataset$tas.retest[i],
        mcq.test = data$dataset$mcq.test[i],
        mcq.retest = data$dataset$mcq.retest[i]
      )

      df_data <- rbind(df_data, row_data)
    }
  }

  return(df_data)
}
