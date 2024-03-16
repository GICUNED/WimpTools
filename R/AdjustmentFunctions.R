
## EMOTIONAL ADJUSTMENT ##

# calculate_adjustment_self_ideal ------------------------------------------------------------
#'
#' This function calculates the Ecuclidean distance and the correlation between the ideal vector and self vector.
#'
#' @param wimp Object containing ideal and self vectors.
#' @param normalize Logical. If TRUE, normalize the distance; otherwise, use the raw Euclidean distance.
#' @return A list with the calculated values: distance and correlation.
#' @export
#'
calculate_adjustment_self_ideal <- function(wimp, normalize = TRUE) {

  # Get ideal and self vectors from the wimp variable
  vector_ideal <- wimp$ideal[[2]]
  vector_self <- wimp$self[[2]]

  # Calculate the euclidean distance between both vectors
  eu_distance <- sqrt(sum((vector_ideal - vector_self)^2))

  # Normalize the Euclidean distance if specified
  distance <- if (normalize) eu_distance / (2 * sqrt(length(vector_ideal)))  else eu_distance

  # Calculate the correlation between the vectors (cosine)
  correlation <- sum(vector_ideal * vector_self) / (sqrt(sum(vector_ideal^2)) * sqrt(sum(vector_self^2)))

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
#' @return The graphic vector for further use.
#'

plot_adjustment_self_ideal <- function(wimp, calculated_values = list()) {

  # Get ideal and self vectors from the wimp variable
  vector_ideal <- wimp$ideal[[2]]
  vector_self <- wimp$self[[2]]

  # Get the calculated graphic vector
  calculated_values <- if (length(calculated_values) == 0) calculate_adjustment_self_ideal(wimp) else calculated_values
  graphic_vector <- c(1, 0, calculated_values$correlation, calculated_values$distance)

  # Set graphic for ideal and self vector
  plot(c(0, vector_ideal[1], vector_self[1]), c(0, vector_ideal[2], vector_self[2]), type = "n", xlab = "Coordenada X", ylab = "Coordenada Y", asp = 1, xlim = c(0, max(vector_ideal[1], vector_self[1]) + 1), ylim = c(0, max(vector_ideal[2], vector_self[2]) + 1))

  arrows(0, 0, vector_ideal[1], vector_ideal[2], col = "green", length = 0.1)
  arrows(0, 0, vector_self[1], vector_self[2], col = "orange", length = 0.1)

  text(vector_ideal[1] + 0.1, vector_ideal[2], "Vector Ideal", col = "green")
  text(vector_self[1] + 0.1, vector_self[2], "Vector Self", col = "orange")

  # Set graphic for adjustment vector: origin (1, 0) and end point (correlation, distance)
  plot(c(0, 1.5), c(0, 1.5), type = "n", xlab = "correlación", ylab = "distancia normalizada", asp = 1, xlim = c(-1, 1), ylim = c(0, 1.5))
  rect(-1, 0, 1, 1, density=20, col='red')
  rect(-0.5, 0, 1, 1, density=20, col='orange')
  rect(0.25, 0, 1, 0.75, density=20, col='yellow')
  rect(0.75, 0, 1, 0.5, density=20, col='green')
  abline(v=c(-1,0,1), lwd = 2, lty=c(2,1,2), col=c("red","black","red"))
  abline(h=c(0,1), lwd = 2, lty=c(1,2), col=c("black","red"))
  arrows(1, 0, calculated_values$correlation, calculated_values$distance, col = "blue", length = 0.1)
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

