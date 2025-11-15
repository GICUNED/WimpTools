## IMPORT FUNCTIONS ##

# Import Weighted Implication Grid ----------------------------------------

#' Import Weighted Implication Grid -- importwimp()
#'
#' @description Function to transform the data of a Weighted Implication
#' Grid (WimpGrid) contained in an Excel file into an S3 object of class
#' wimp. This function reads the Excel data, normalizes values, calculates
#' hypothetical scenarios, and structures the data into vertices and edges
#' for network analysis.
#'
#' @param path Path to the Excel file on your computer. The file suffix
#'        must be .xlsx.
#' @param sheet Number of the Excel sheet that contains the WimpGrid data.
#'        Defaults to 1.
#'
#' @return A wimp S3 object containing:
#'   \itemize{
#'     \item global: metadata and matrices (scale, wmatrix, hypo_matrix)
#'     \item vertices: data.frame with construct information
#'     \item edges: data.frame with implication relationships
#'   }
#'
#' @examples
#' \dontrun{
#' # Import your own Excel file
#' my_wimp <- importwimp("path/to/your/file.xlsx")
#' # Or use the example data included in the package
#' example_path <- system.file("extdata", "example_wimp.xlsx",
#'                            package = "WimpTools")
#' example_wimp <- importwimp(example_path)
#' #' }
#'
#' @export
#'
#' @import readxl

importwimp <- function(path, sheet = 1) {

  # Read Excel file and suppress column name messages
  xlsx <- suppressMessages(readxl::read_excel(path, sheet = sheet,
                                              col_names = FALSE))
  global_list <- list()

  # Extract global metadata from columns 1-2 (starting from row 2)
  for (r in 2:nrow(xlsx)){
    col1_val <- xlsx[[r, 1]]
    col2_val <- xlsx[[r, 2]]
    if (!is.na(col1_val) && as.character(col1_val) != "") {
      val <- if (!is.na(col2_val) && as.character(col2_val) != "") {
        as.character(col2_val)
      } else {
        ""
      }
      global_list[[as.character(col1_val)]] <- val
    }
  }

  # Determine number of constructs from Excel structure
  n_constructs <- nrow(xlsx[3]) - 1
  if (n_constructs < 1) stop("No constructs found in Excel file.")
  row_constructs <- 2:(n_constructs + 1)

  # Extract data from Excel columns:
  # Column 3: left poles, Column 4: self ratings,
  # Columns 5 to (4+n): hypothetical matrix,
  # Column (5+n): ideal ratings, Column (6+n): right poles
  left_poles <- as.character(xlsx[row_constructs, 3][[1]])
  self_col <- 4
  self_vector <- as.numeric(xlsx[row_constructs, self_col][[1]])
  # Extract hypothetical matrix (n x n square matrix)
  # Rows represent constructs, columns represent hypothetical scenarios
  hypo_start_col <- self_col + 1
  hypo_end_col <- hypo_start_col + n_constructs - 1
  hypo_matrix_raw <- xlsx[row_constructs, hypo_start_col:hypo_end_col]
  hypo_matrix <- as.matrix(sapply(hypo_matrix_raw, as.numeric))
  diag(hypo_matrix) <- self_vector  # Replace diagonal with actual self ratings
  # Extract scale range and ideal ratings
  scale_min <- as.numeric(xlsx[[1, 3]])
  ideal_col_idx <- hypo_end_col + 1
  scale_max_col <- ideal_col_idx + 1
  scale_max <- as.numeric(xlsx[[1, scale_max_col]])
  scale_center <- (scale_min + scale_max) / 2
  direct_ideal <- as.numeric(xlsx[row_constructs, ideal_col_idx][[1]])
  right_poles <- as.character(xlsx[row_constructs, scale_max_col][[1]])
  # Extra attributes (if any)
  extra_attrs_start <- scale_max_col + 1
  extra_attrs_end <- ncol(xlsx)
  extra_attrs <- if (extra_attrs_start <= extra_attrs_end) {
    xlsx[row_constructs, extra_attrs_start:extra_attrs_end]
  } else {
    NULL
  }
  extra_attr_names <- if (!is.null(extra_attrs) && ncol(extra_attrs) > 0) {
    unlist(as.list(xlsx[1, extra_attrs_start:extra_attrs_end]))
  } else {
    NULL
  }

  # Normalize all ratings to [-1, 1] range using scale center
  normalized_self <- (self_vector - (scale_center * rep(1, n_constructs))) /
    (0.5 * (scale_max - scale_min))
  normalized_ideal <- (direct_ideal - (scale_center * rep(1, n_constructs))) /
    (0.5 * (scale_max - scale_min))
  # Calculate hypothetical values using .calc.hypo algorithm
  normalized_hypothetical <- mapply(.calc.hypo, normalized_self,
                                    normalized_ideal)
  # Normalize the entire hypothetical matrix and update diagonal
  normalized_hypo_matrix <- (hypo_matrix -
                               (scale_center *
                                  matrix(rep(1, n_constructs * n_constructs),
                                         ncol = n_constructs))) /
    (0.5 * (scale_max - scale_min))
  diag(normalized_hypo_matrix) <- normalized_hypothetical
  # Generate descriptive names for hypothetical scenarios
  # Names follow pattern "Totally [Pole]" based on hypothetical direction
  hypo_names <- character(n_constructs)
  for (i in seq_len(n_constructs)) {
    hypo_val <- normalized_hypothetical[i]
    if (!is.na(hypo_val)) {
      if (hypo_val > 0) {
        hypo_names[i] <- paste("Totally", right_poles[i])
      } else {
        hypo_names[i] <- paste("Totally", left_poles[i])
      }
    } else {
      hypo_names[i] <- paste("Totally", left_poles[i], "-", right_poles[i])
    }
  }
  colnames(normalized_hypo_matrix) <- hypo_names
  rownames(normalized_hypo_matrix) <- paste(left_poles, "-", right_poles)
  self_poles <- mapply(.self.poles, normalized_self, left_poles, right_poles)
  ideal_poles <- mapply(.self.poles, normalized_ideal, left_poles, right_poles)

  # Calculate weight matrix for implication analysis
  # This represents the standardized implication relationships
  imp_matrix_norm <- (hypo_matrix -
                        (scale_center *
                           matrix(rep(1, n_constructs * n_constructs),
                                  ncol = n_constructs))) /
    (0.5 * (scale_max - scale_min))
  imp_matrix <- t(imp_matrix_norm)
  num_weight_matrix <- imp_matrix - matrix(normalized_self,
                                           nrow = n_constructs,
                                           ncol = n_constructs,
                                           byrow = TRUE)
  den_weight_matrix <- matrix(normalized_hypothetical,
                              nrow = n_constructs,
                              ncol = n_constructs) -
    matrix(normalized_self,
           nrow = n_constructs,
           ncol = n_constructs)
  weight_matrix <- num_weight_matrix / den_weight_matrix
  # Assign construct names to weight matrix rows and columns
  construct_names <- paste(left_poles, "-", right_poles)
  colnames(weight_matrix) <- construct_names
  rownames(weight_matrix) <- construct_names
  # Construct the wimp S3 object
  wimp <- list()
  class(wimp) <- c("wimp", "list")

  # Store global metadata and matrices
  global_list$scale <- c(scale_min, scale_max)
  global_list$n_constructs <- as.numeric(n_constructs)
  global_list$weight_matrix <- weight_matrix
  global_list$hypo_matrix <- normalized_hypo_matrix

  # Classify construct congruence types based on self-ideal relationships
  congruence <- character(n_constructs)
  for (k in seq_len(n_constructs)) {
    if (is.na(normalized_ideal[k]) || normalized_ideal[k] == 0) {
      congruence[k] <- "Dilemmatic"
    } else if (is.na(normalized_self[k]) || normalized_self[k] == 0) {
      congruence[k] <- "Undefined"
    } else if (sign(normalized_self[k]) == sign(normalized_ideal[k])) {
      congruence[k] <- "Congruent"
    } else {
      congruence[k] <- "Discrepant"
    }
  }
  # Create vertices data.frame: one row per construct with all attributes
  vertices_df <- data.frame(
    id = seq_len(n_constructs),
    left_pole = left_poles,
    right_pole = right_poles,
    self = normalized_self,
    ideal = normalized_ideal,
    self_pole = self_poles,
    ideal_pole = ideal_poles,
    congruence = congruence,
    stringsAsFactors = FALSE
  )
  # Add extra attributes from Excel if present
  if (!is.null(extra_attrs) && ncol(extra_attrs) > 0) {
    for (col_idx in seq_len(ncol(extra_attrs))) {
      col_name <- if (!is.null(extra_attr_names) &&
                        col_idx <= length(extra_attr_names)) {
        extra_attr_names[col_idx]
      } else {
        paste0("attr_", col_idx)
      }
      if (is.na(col_name) || col_name == "") {
        col_name <- paste0("attr_", col_idx)
      }
      vertices_df[[col_name]] <- as.character(extra_attrs[[col_idx]])
    }
  }
  # Create edges data.frame: one row per significant implication weight
  # Excludes self-loops and near-zero weights
  edges_list <- vector("list", 0)
  for (i in seq_len(n_constructs)) {
    for (j in seq_len(n_constructs)) {
      if (i == j) next  # Skip self-loops
      wval <- as.numeric(weight_matrix[i, j])
      if (is.na(wval) || abs(wval) < .Machine$double.eps) next
      edges_list[[length(edges_list) + 1]] <- list(
        id = paste(i, "t", j, sep = ""),
        from = as.integer(i),
        to = as.integer(j),
        weight = wval
      )
    }
  }

  # Convert edge list to data.frame format

  edges_df <- if (length(edges_list) > 0) {
    tmp <- do.call(rbind, lapply(edges_list,
                                 function(z) {
                                   as.data.frame(z, stringsAsFactors = FALSE)
                                 }))
    tmp$from <- as.integer(tmp$from)
    tmp$to <- as.integer(tmp$to)
    tmp$weight <- as.numeric(tmp$weight)
    tmp
  } else {
    data.frame(id = character(), from = integer(), to = integer(),
               weight = numeric(), stringsAsFactors = FALSE)
  }
  # Assemble final wimp object with all components
  wimp$global <- global_list
  wimp$vertices <- vertices_df
  wimp$edges <- edges_df
  return(wimp)
}