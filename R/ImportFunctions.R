## IMPORT FUNCTIONS ##

# Import Weigthed ImpGrid -------------------------------------------------

#' Import Weighted Implication Grid -- importwimp()
#'
#' @description Function to transform the data of a WimpGrid contained in an
#' Excel file into an S3 object of class wimp.
#'
#' @param path Path to the excel file on your computer. The file suffix
#'        has to be .xlsx.
#' @param sheet Number of the Excel sheet that contains the WimpGrid data.
#'
#' @return A wimp S3 object.
#'
#' @export
#'
#' @import readxl

importwimp <- function(path, sheet = 1) {

  # Read Excel file
  xlsx <- suppressMessages(readxl::read_excel(path, sheet = sheet,
                                              col_names = FALSE))
  global_list <- list()

  # Extract global attributes from columns 1-2 (rows 2+)
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

  # Number of constructs = rows - 1 (row 1 is header)
  n_constructs <- nrow(xlsx[3]) - 1
  if (n_constructs < 1) stop("No constructs found in Excel file.")
  row_constructs <- 2:(n_constructs + 1)

  # Column layout: 3=lpole, 4=self, 5:(4+n)=hypo,
  # (5+n)=ideal, (6+n)=rpole, (7+n)+extras
  left_poles <- as.character(xlsx[row_constructs, 3][[1]])
  self_col <- 4
  self_vector <- as.numeric(xlsx[row_constructs, self_col][[1]])
  # Hypothetical matrix (n x n, rows=constructs, cols=constructs)
  hypo_start_col <- self_col + 1
  hypo_end_col <- hypo_start_col + n_constructs - 1
  hypo_matrix_raw <- xlsx[row_constructs, hypo_start_col:hypo_end_col]
  hypo_matrix <- as.matrix(sapply(hypo_matrix_raw, as.numeric))
  diag(hypo_matrix) <- self_vector  # Replace header labels with actual values
  # Scale and ideal values
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

  # Normalize to (-1, 1) range
  normalized_self <- (self_vector - (scale_center * rep(1, n_constructs))) /
    (0.5 * (scale_max - scale_min))
  normalized_ideal <- (direct_ideal - (scale_center * rep(1, n_constructs))) /
    (0.5 * (scale_max - scale_min))
  normalized_hypothetical <- mapply(.calc.hypo, normalized_self,
                                    normalized_ideal)
  self_poles <- mapply(.self.poles, normalized_self, left_poles, right_poles)
  ideal_poles <- mapply(.self.poles, normalized_ideal, left_poles, right_poles)

  # Calculate weight matrix (standardized implications)
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
  # Build wimp object
  wimp <- list()
  class(wimp) <- c("wimp", "list")

  # Global: metadata
  global_list$scale <- c(scale_min, scale_max)
  global_list$n_constructs <- as.numeric(n_constructs)
  global_list$wmatrix <- weight_matrix

  # Congruency classification
  congruency <- character(n_constructs)
  for (k in seq_len(n_constructs)) {
    if (is.na(normalized_ideal[k]) || normalized_ideal[k] == 0) {
      congruency[k] <- "dilemmatic"
    } else if (is.na(normalized_self[k]) || normalized_self[k] == 0) {
      congruency[k] <- "undefined"
    } else if (sign(normalized_self[k]) == sign(normalized_ideal[k])) {
      congruency[k] <- "congruent"
    } else {
      congruency[k] <- "discrepant"
    }
  }
  # Vertices: one row per construct
  vertices_df <- data.frame(
    id = seq_len(n_constructs),
    lpole = left_poles,
    rpole = right_poles,
    self = normalized_self,
    ideal = normalized_ideal,
    self_pole = self_poles,
    ideal_pole = ideal_poles,
    congruency = congruency,
    stringsAsFactors = FALSE
  )
  # Append extra attributes to vertices
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
  # Edges: one row per non-zero weight (excluding self-loops)
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
  wimp$global <- global_list
  wimp$vertices <- vertices_df
  wimp$edges <- edges_df
  return(wimp)
}
