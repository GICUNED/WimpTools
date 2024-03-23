
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

calculate_adjustment_self_ideal <- function(wimp, normalize = TRUE, filtered_constructs) {

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
  magnitude <- if (normalize) 1 - (sqrt(correlation^2 + distance^2)) / sqrt(5) else NULL

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

plot_adjustment_self_ideal <- function(wimp, calculated_values = list(), filtered_constructs = list()) {

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

Adjustment_self_ideal <- function(wimp, normalize = TRUE) {
  # Calculate adjustment and self-ideal values
  calculated_values <- calculate_adjustment_self_ideal(wimp, normalize = normalize)

  # Plot the graphs
  graphic_vector <- plot_adjustment_self_ideal(wimp, calculated_values)

  # Return data
  return(list(distance = calculated_values$distance, correlation = calculated_values$correlation, magnitude = calculated_values$magnitude, graphic = graphic_vector))
}

