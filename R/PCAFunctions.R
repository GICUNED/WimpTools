## WimpGrid Biplot Function

#' Weighted Implication Grid Biplot — wimp_biplot()
#'
#' @description PCA biplot combining hypothetical selves and personal
#' constructs. Visualizes the first two principal components to reveal
#' patterns and relationships among constructs and situations.
#'
#' @param wimp A `wimp` object imported via \code{importwimp()}.
#' @param text_size Text label size multiplier. Defaults to 1.
#'
#' @return A plotly biplot with PC1/PC2 axes, construct vectors, and
#'   situation labels (SELF, IDEAL, hypothetical scenarios).
#'
#' @details Constructs are reoriented so all ideal values are positive
#'   before PCA. Construct poles are displayed as vectors radiating from
#'   the origin.
#'
#' @author Alejandro Sanfeliciano
#' @importFrom DescTools CartToPol PolToCart
#' @import useful
#' @importFrom stats cor prcomp runif
#' @export
#' @examples
#' wimp_biplot(example_wimp)

wimp_biplot <- function(wimp, text_size = 1) {
  # Align constructs so all ideal values are positive
  wimp <- .align_wimp(wimp, exclude_dilemmatics = FALSE)

  left_poles <- wimp$vertices$left_pole
  right_poles <- wimp$vertices$right_pole
  n_constructs <- length(left_poles)

  # Build matrix: SELF | hypo scenarios | IDEAL
  hypo_matrix <- wimp$global$hypo_matrix
  self_vector <- wimp$vertices$self
  ideal_vector <- wimp$vertices$ideal
  pca_matrix <- cbind(self_vector, hypo_matrix, ideal_vector)
  colnames(pca_matrix)[c(1, ncol(pca_matrix))] <- c("SELF", "IDEAL")

  # PCA on first two principal components
  pca <- prcomp(pca_matrix, rank. = 2)

  # Variance explained percentages
  variance_pc1 <- round(
    summary(pca)[6]$importance[2, 1] * 100, digits = 2
  )
  variance_pc2 <- round(
    summary(pca)[6]$importance[2, 2] * 100, digits = 2
  )

  # Rotation coordinates (construct loadings)
  pc1_loads <- pca$rotation[, 1]
  pc2_loads <- pca$rotation[, 2]

  # Convert to polar for vector placement
  pc_radius <- DescTools::CartToPol(pc1_loads, pc2_loads)$r

  situation_names <- colnames(pca_matrix)

  # Construct vector placement (opposite poles)
  construct_x <- c(pca$x[, 1], -pca$x[, 1])
  construct_y <- c(pca$x[, 2], -pca$x[, 2])

  # Add spacing to avoid overlap
  vector_radius <- max(abs(pc_radius)) + 0.05 * max(abs(pc_radius)) +
    max(abs(pc_radius)) * runif(n_constructs, -0.35, 0.1)
  vector_theta <- DescTools::CartToPol(construct_x, construct_y)$theta

  vector_x <- DescTools::PolToCart(vector_radius, vector_theta)$x
  vector_y <- DescTools::PolToCart(vector_radius, vector_theta)$y

  pole_names <- c(right_poles, left_poles)

  # Prepare data frames for plotting
  df_situations <- data.frame(situation_names, pc1_loads, pc2_loads)
  df_vectors <- data.frame(pole_names, vector_x, vector_y)

  # Plot range
  plot_range <- max(abs(df_vectors[2:3]) + 0.1 * abs(df_vectors[2:3]))

  # Build plotly biplot
  fig <- plot_ly(type = "scatter") %>%
    # Hypothetical situations (middle rows, excluding SELF and IDEAL)
    add_annotations(
      data = df_situations[-c(1, n_constructs + 2), ],
      x = ~pc1_loads, y = ~pc2_loads, text = ~situation_names,
      hoverinfo = "text", font = list(size = 15 * text_size),
      showarrow = FALSE, xanchor = "center", yanchor = "center"
    ) %>%
    # SELF (first row, dark blue)
    add_annotations(
      data = df_situations[1, ],
      x = ~pc1_loads, y = ~pc2_loads, text = ~situation_names,
      hoverinfo = "text",
      font = list(size = 15 * text_size, color = "darkblue"),
      showarrow = FALSE, xanchor = "center", yanchor = "center"
    ) %>%
    # IDEAL (last row, dark green)
    add_annotations(
      data = df_situations[n_constructs + 2, ],
      x = ~pc1_loads, y = ~pc2_loads, text = ~situation_names,
      hoverinfo = "text",
      font = list(size = 15 * text_size, color = "darkgreen"),
      showarrow = FALSE, xanchor = "center", yanchor = "center"
    ) %>%
    # Construct vectors as dashed lines
    add_segments(
      data = df_vectors, x = 0, xend = ~vector_x, y = 0, yend = ~vector_y,
      line = list(color = "#6F6BFF", dash = "dot", width = 0.75),
      hoverinfo = "none", inherit = FALSE, showlegend = FALSE
    ) %>%
    # Construct pole labels
    add_annotations(
      data = df_vectors, x = ~vector_x, y = ~vector_y, text = ~pole_names,
      hoverinfo = "none",
      font = list(size = 12 * text_size, color = "#6F6BFF"),
      showarrow = FALSE, xanchor = "center", yanchor = "bottom"
    ) %>%
    layout(
      xaxis = list(
        title = paste("<B>PC1</B> [", variance_pc1, "%]", sep = ""),
        range = c(-plot_range, plot_range), zeroline = TRUE,
        zerolinecolor = "black", zerolinewidth = 2, showline = FALSE
      ),
      yaxis = list(
        title = paste("<B>PC2</B> [", variance_pc2, "%]", sep = ""),
        range = c(-plot_range, plot_range), zeroline = TRUE,
        zerolinecolor = "black", zerolinewidth = 2, showline = FALSE
      ),
      showlegend = FALSE
    )

  return(fig)
}
