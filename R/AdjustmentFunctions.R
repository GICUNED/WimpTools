## EMOTIONAL ADJUSTMENT ##

#' Title
#'
#' @param wimp
#'
#' @return
#' @export
#'
#' @examples
Adjustment_self_ideal <- function(wimp) {
  # Get ideal and self vectors from wimp variable
  vector_ideal <- wimp$ideal[[2]]
  vector_self <- wimp$self[[2]]

  # Calculate the euclidean distance between both vectors
  eu_distance <- sqrt(sum((vector_ideal - vector_self)^2))

  # Normalize euclidean distance
  normalized_distance <- eu_distance / (2 * sqrt(length(vector_ideal)))

  # Calculate the correlation between the vectors (cosine)
  correlation <- sum(vector_ideal * vector_self) / (sqrt(sum(vector_ideal^2)) * sqrt(sum(vector_self^2)))

  #Return data
  return(list(normalized_distance = normalized_distance, correlation = correlation))
}


#' Title
#'
#' @return
#' @export
#'
#' @examples
graph_adjustment <- function(wimp) {

  graphic_vector_data <- Adjustment_self_ideal(wimp)

  # Set up the graphic vector with origin (1,0)
  graphic_vector <- c(1, 0, graphic_vector_data$correlation, graphic_vector_data$normalized_distance)

  # Set graphic for ideal and self vectors
  plot(c(0, vector_ideal[1], vector_self[1]), c(0, vector_ideal[2], vector_self[2]), type = "n", xlab = "Coordenada X", ylab = "Coordenada Y", asp = 1, xlim = c(0, max(vector_ideal[1], vector_self[1]) + 1), ylim = c(0, max(vector_ideal[2], vector_self[2]) + 1))
  arrows(0, 0, vector_ideal[1], vector_ideal[2], col = "green", length = 0.1)
  arrows(0, 0, vector_self[1], vector_self[2], col = "orange", length = 0.1)
  text(vector_ideal[1] + 0.1, vector_ideal[2], "Vector Ideal", col = "green")
  text(vector_self[1] + 0.1, vector_self[2], "Vector Self", col = "orange")

  # Set graphic for adjustment vector: origin (1, 0) and end point (correlation, normalized_distance)
  plot(c(0, 1.5), c(0, 1.5), type = "n", xlab = "Correlation", ylab = "Normalized Distance", asp = 1, xlim = c(0, 1.5), ylim = c(0, 1.5))
  arrows(1, 0, correlation, normalized_distance, col = "blue", length = 0.1)
  text(1.2, 0.1, "Vector (1,0) to (correlation, normalized_distance)", col = "blue")
}
