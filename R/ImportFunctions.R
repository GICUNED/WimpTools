## IMPORT FUNCTIONS ##

utils::globalVariables(c(".self_poles"))

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
#' @references
#' Sanfeliciano, A., et al. (2024). The Weighted Implications Grid: A Graph-Theoretical Approach to Modelling Psychological Change Construction.
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
  normalized_hypothetical <- mapply(.calc_hypo, normalized_self,
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
  self_poles <- mapply(.self_poles, normalized_self, left_poles, right_poles)
  ideal_poles <- mapply(.self_poles, normalized_ideal, left_poles, right_poles)

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
      
      # Special handling for 'preference' column: normalize using wimp scale
      if (tolower(col_name) == "preference") {
        pref_raw <- as.numeric(extra_attrs[[col_idx]])
        # Normalize to [-1, 1] using the same scale as self and ideal
        pref_normalized <- (pref_raw - (scale_center * rep(1, n_constructs))) /
          (0.5 * (scale_max - scale_min))
        vertices_df[[col_name]] <- pref_normalized
      } else {
        vertices_df[[col_name]] <- as.character(extra_attrs[[col_idx]])
      }
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


# Import WimpGrid from a JSON string or R list -----------------------------

#' Import WimpGrid from a JSON string or R list
#'
#' @description Function to transform the data of a Weighted Implication
#' Grid (WimpGrid) in JSON format (or pre-parsed R list) into an S3 object
#' of class wimp. This supports integration with database web app workers,
#' such as psychlab.
#'
#' @param json_data A JSON character string, a file path to a JSON file, or
#'        an already parsed list representing the grid data.
#'
#' @return A wimp S3 object containing:
#'   \itemize{
#'     \item global: metadata and matrices (scale, weight_matrix, hypo_matrix)
#'     \item vertices: data.frame with construct information
#'     \item edges: data.frame with implication relationships
#'   }
#'
#' @export
#'
#' @importFrom jsonlite fromJSON
importwimp_json <- function(json_data) {

  # 1. Parse JSON if string or load list
  if (is.character(json_data)) {
    if (length(json_data) == 1 && file.exists(json_data)) {
      parsed <- jsonlite::fromJSON(json_data, simplifyVector = TRUE)
    } else {
      parsed <- jsonlite::fromJSON(json_data, simplifyVector = TRUE)
    }
  } else if (is.list(json_data)) {
    parsed <- json_data
  } else {
    stop("Input must be a JSON character string, a file path to a JSON file, or a parsed list.")
  }

  # Helper for robust name matching
  get_field <- function(df_or_list, possible_names, required = FALSE) {
    names <- names(df_or_list)
    found_name <- NULL
    for (name in possible_names) {
      match_idx <- which(tolower(names) == tolower(name))
      if (length(match_idx) > 0) {
        found_name <- names[match_idx[1]]
        break
      }
    }
    if (is.null(found_name)) {
      if (required) {
        stop(paste("Required field missing. Tried matching names:", paste(possible_names, collapse = ", ")))
      }
      return(NULL)
    }
    return(df_or_list[[found_name]])
  }

  # 2. Extract and Validate Scale
  scale_val <- get_field(parsed, c("scale", "range"), required = TRUE)
  if (is.list(scale_val) && !is.numeric(scale_val)) {
    scale_min <- as.numeric(get_field(scale_val, c("min", "scale_min", "minimum"), required = TRUE))
    scale_max <- as.numeric(get_field(scale_val, c("max", "scale_max", "maximum"), required = TRUE))
  } else if (is.numeric(scale_val) && length(scale_val) >= 2) {
    scale_min <- as.numeric(scale_val[1])
    scale_max <- as.numeric(scale_val[2])
  } else {
    stop("scale must be an object with 'min' and 'max' fields, or a numeric vector of length 2.")
  }

  # 3. Extract and Parse Constructs
  constructs_val <- get_field(parsed, c("constructs", "nodes", "variables"), required = TRUE)
  
  # Pre-process 'attrs' if they exist (flatten nested attributes)
  if (is.data.frame(constructs_val)) {
    if ("attrs" %in% names(constructs_val) && is.data.frame(constructs_val$attrs)) {
      attrs_df <- constructs_val$attrs
      constructs_val$attrs <- NULL
      constructs_val <- cbind(constructs_val, attrs_df)
    }
  } else if (is.list(constructs_val)) {
    constructs_val <- lapply(constructs_val, function(item) {
      if ("attrs" %in% names(item) && is.list(item[["attrs"]])) {
        attrs_list <- item[["attrs"]]
        item[["attrs"]] <- NULL
        item <- c(item, attrs_list)
      }
      return(item)
    })
  }
  
  if (is.list(constructs_val) && !is.data.frame(constructs_val)) {
    # Nested list format (array of objects)
    n_constructs <- length(constructs_val)
    if (n_constructs < 1) stop("No constructs found in constructs list.")
    
    left_poles <- character(n_constructs)
    right_poles <- character(n_constructs)
    self_vector <- numeric(n_constructs)
    direct_ideal <- numeric(n_constructs)
    
    all_keys <- unique(unlist(lapply(constructs_val, names)))
    known_keys <- c("left_pole", "leftpole", "left_poles", "lpole", "left",
                    "right_pole", "rightpole", "right_poles", "rpole", "right",
                    "self", "self_rating", "self_vector",
                    "ideal", "ideal_rating", "ideal_vector", "direct_ideal")
    extra_keys <- setdiff(all_keys, known_keys)
    
    extra_attrs <- list()
    for (k in extra_keys) {
      extra_attrs[[k]] <- vector("character", n_constructs)
    }
    
    for (i in seq_len(n_constructs)) {
      item <- constructs_val[[i]]
      left_poles[i] <- as.character(get_field(item, c("left_pole", "leftpole", "left_poles", "lpole", "left"), required = TRUE))
      right_poles[i] <- as.character(get_field(item, c("right_pole", "rightpole", "right_poles", "rpole", "right"), required = TRUE))
      self_vector[i] <- as.numeric(get_field(item, c("self", "self_rating", "self_vector"), required = TRUE))
      direct_ideal[i] <- as.numeric(get_field(item, c("ideal", "ideal_rating", "ideal_vector", "direct_ideal"), required = TRUE))
      
      for (k in extra_keys) {
        val <- item[[k]]
        extra_attrs[[k]][i] <- if (is.null(val) || is.na(val)) NA_character_ else as.character(val)
      }
    }
    
    constructs_df <- data.frame(
      left_pole = left_poles,
      right_pole = right_poles,
      self = self_vector,
      ideal = direct_ideal,
      stringsAsFactors = FALSE
    )
    if (length(extra_keys) > 0) {
      for (k in extra_keys) {
        constructs_df[[k]] <- extra_attrs[[k]]
      }
    }
  } else if (is.data.frame(constructs_val)) {
    # Simplified data frame format
    constructs_df <- constructs_val
    n_constructs <- nrow(constructs_df)
    if (n_constructs < 1) stop("No constructs found in constructs data.frame.")
    
    left_poles <- as.character(get_field(constructs_df, c("left_pole", "leftpole", "left_poles", "lpole", "left"), required = TRUE))
    right_poles <- as.character(get_field(constructs_df, c("right_pole", "rightpole", "right_poles", "rpole", "right"), required = TRUE))
    self_vector <- as.numeric(get_field(constructs_df, c("self", "self_rating", "self_vector"), required = TRUE))
    direct_ideal <- as.numeric(get_field(constructs_df, c("ideal", "ideal_rating", "ideal_vector", "direct_ideal"), required = TRUE))
  } else {
    stop("constructs must be a list of construct objects or a data.frame.")
  }

  # Range Validation
  if (any(self_vector < scale_min | self_vector > scale_max, na.rm = TRUE)) {
    warning("Some 'self' ratings are out of scale bounds.")
  }
  if (any(direct_ideal < scale_min | direct_ideal > scale_max, na.rm = TRUE)) {
    warning("Some 'ideal' ratings are out of scale bounds.")
  }

  # 4. Extract and Validate Hypothetical Matrix
  raw_hypo <- get_field(parsed, c("hypo_matrix", "hypo_grid", "hypomatrix", "hypo", "imp_matrix"), required = TRUE)
  if (is.list(raw_hypo) && !is.matrix(raw_hypo)) {
    hypo_matrix <- do.call(rbind, lapply(raw_hypo, as.numeric))
  } else if (is.matrix(raw_hypo)) {
    hypo_matrix <- apply(raw_hypo, 2, as.numeric)
  } else {
    stop("hypo_matrix must be a matrix or a list of rows.")
  }
  
  if (nrow(hypo_matrix) != n_constructs || ncol(hypo_matrix) != n_constructs) {
    stop(sprintf("hypo_matrix must be a square matrix of size %d x %d. Got %d x %d.", 
                 n_constructs, n_constructs, nrow(hypo_matrix), ncol(hypo_matrix)))
  }
  
  diag(hypo_matrix) <- self_vector  # Replace diagonal with actual self ratings

  # 5. Global Metadata Extraction
  global_list <- list()
  meta <- get_field(parsed, c("metadata", "meta", "global_attrs"))
  if (is.list(meta)) {
    for (name in names(meta)) {
      val <- meta[[name]]
      global_list[[name]] <- if (is.null(val) || is.na(val)) "" else as.character(val)
    }
  }

  # 6. Normalize Ratings
  scale_center <- (scale_min + scale_max) / 2
  normalized_self <- (self_vector - (scale_center * rep(1, n_constructs))) /
    (0.5 * (scale_max - scale_min))
  normalized_ideal <- (direct_ideal - (scale_center * rep(1, n_constructs))) /
    (0.5 * (scale_max - scale_min))
  
  normalized_hypothetical <- mapply(.calc_hypo, normalized_self, normalized_ideal)
  
  normalized_hypo_matrix <- (hypo_matrix -
                               (scale_center *
                                  matrix(rep(1, n_constructs * n_constructs),
                                         ncol = n_constructs))) /
    (0.5 * (scale_max - scale_min))
  diag(normalized_hypo_matrix) <- normalized_hypothetical
  
  # Generate descriptive scenario names
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
  
  self_poles <- mapply(.self_poles, normalized_self, left_poles, right_poles)
  ideal_poles <- mapply(.self_poles, normalized_ideal, left_poles, right_poles)

  # 7. Calculate weight matrix for implication analysis
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
  
  construct_names <- paste(left_poles, "-", right_poles)
  colnames(weight_matrix) <- construct_names
  rownames(weight_matrix) <- construct_names

  # 8. Construct Classifications
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

  # 9. Create vertices data.frame
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
  
  # Populate any extra columns
  known_cols <- c("left_pole", "leftpole", "left_poles", "lpole", "left",
                  "right_pole", "rightpole", "right_poles", "rpole", "right",
                  "self", "self_rating", "self_vector",
                  "ideal", "ideal_rating", "ideal_vector", "direct_ideal")
  extra_cols <- setdiff(colnames(constructs_df), known_cols)
  if (length(extra_cols) > 0) {
    for (col_name in extra_cols) {
      if (tolower(col_name) == "preference") {
        pref_raw <- as.numeric(constructs_df[[col_name]])
        pref_normalized <- (pref_raw - (scale_center * rep(1, n_constructs))) /
          (0.5 * (scale_max - scale_min))
        vertices_df[[col_name]] <- pref_normalized
      } else {
        vertices_df[[col_name]] <- as.character(constructs_df[[col_name]])
      }
    }
  }

  # 10. Create edges list
  edges_list <- vector("list", 0)
  for (i in seq_len(n_constructs)) {
    for (j in seq_len(n_constructs)) {
      if (i == j) next
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

  # 11. Assemble S3 object
  wimp <- list()
  class(wimp) <- c("wimp", "list")
  
  global_list$scale <- c(scale_min, scale_max)
  global_list$n_constructs <- as.numeric(n_constructs)
  global_list$weight_matrix <- weight_matrix
  global_list$hypo_matrix <- normalized_hypo_matrix
  
  wimp$global <- global_list
  wimp$vertices <- vertices_df
  wimp$edges <- edges_df
  return(wimp)
}