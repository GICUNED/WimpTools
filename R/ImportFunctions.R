## IMPORT FUNCTIONS ##

# Import Weigthed ImpGrid -------------------------------------------------

#' Import Weighted Implication Grid -- importwimp()
#'
#' @description Function to transform the data of a WimpGrid contained in an
#' Excel file into an S3 object of class wimp.
#'
#' @param path Path to the excel file on your computer. The file suffix has to be .xlsx.
#' @param sheet Number of the Excel sheet that contains the WimpGrid data.
#'
#' @return A wimp S3 object.
#'
#' @export
#'
#' @import readxl


importwimp <- function(path, sheet = 1){

  # Read Excel file
  xlsx <- suppressMessages(readxl::read_excel(path, sheet = sheet, col_names = FALSE))
  
  global_list <- list()

  # Extract global attributes from columns 1-2 (rows 2+)
  for(r in 2:nrow(xlsx)){
    col1_val <- xlsx[[r, 1]]
    col2_val <- xlsx[[r, 2]]
    if(!is.na(col1_val) && as.character(col1_val) != ""){
      global_list[[ as.character(col1_val) ]] <- if(!is.na(col2_val) && as.character(col2_val) != "") as.character(col2_val) else ""
    }
  }

  # Number of constructs = rows - 1 (row 1 is header)
  n.constructs <- nrow(xlsx[3]) - 1
  if(n.constructs < 1) stop("No constructs found in Excel file.")
  row.constructs <- 2:(n.constructs + 1)

  # Column layout: 3=lpole, 4=self, 5:(4+n)=hypo, (5+n)=ideal, (6+n)=rpole, (7+n)+extras
  left.poles <- as.character(xlsx[row.constructs, 3][[1]])
  self.col <- 4
  self.vector <- as.numeric(xlsx[row.constructs, self.col][[1]])
  
  # Hypothetical matrix (n x n, rows=constructs, cols=constructs)
  hypo_start_col <- self.col + 1
  hypo_end_col <- hypo_start_col + n.constructs - 1
  hypo_matrix_raw <- xlsx[row.constructs, hypo_start_col:hypo_end_col]
  hypo_matrix <- as.matrix(sapply(hypo_matrix_raw, as.numeric))
  diag(hypo_matrix) <- self.vector  # Replace header labels with actual self values
  
  # Scale and ideal values
  scale.min <- as.numeric(xlsx[[1, 3]])
  ideal_col_idx <- hypo_end_col + 1
  scale.max_col <- ideal_col_idx + 1
  scale.max <- as.numeric(xlsx[[1, scale.max_col]])
  scale.center <- (scale.min + scale.max) / 2
  direct.ideal <- as.numeric(xlsx[row.constructs, ideal_col_idx][[1]])
  right.poles <- as.character(xlsx[row.constructs, scale.max_col][[1]])
  
  # Extra attributes (if any)
  extra_attrs_start <- scale.max_col + 1
  extra_attrs_end <- ncol(xlsx)
  extra_attrs <- if(extra_attrs_start <= extra_attrs_end) xlsx[row.constructs, extra_attrs_start:extra_attrs_end] else NULL
  extra_attr_names <- if(!is.null(extra_attrs) && ncol(extra_attrs) > 0) 
                        unlist(as.list(xlsx[1, extra_attrs_start:extra_attrs_end])) else NULL
  
  # Normalize to (-1, 1) range
  normalized.self <- (self.vector - (scale.center * rep(1, n.constructs))) / (0.5 * (scale.max - scale.min))
  normalized.ideal <- (direct.ideal - (scale.center * rep(1, n.constructs))) / (0.5 * (scale.max - scale.min))
  normalized.hypothetical <- mapply(.calc.hypo, normalized.self, normalized.ideal)
  self.poles <- mapply(.self.poles, normalized.self, left.poles, right.poles)
  ideal.poles <- mapply(.self.poles, normalized.ideal, left.poles, right.poles)
  
  # Calculate weight matrix (standardized implications)
  imp.matrix_norm <- (hypo_matrix - (scale.center * matrix(rep(1, n.constructs * n.constructs), ncol = n.constructs))) / 
                     (0.5 * (scale.max - scale.min))
  imp.matrix <- t(imp.matrix_norm)
  num.weight.matrix <- imp.matrix - matrix(normalized.self, nrow = n.constructs, ncol = n.constructs, byrow = TRUE)
  den.weigth.matrix <- matrix(normalized.hypothetical, nrow = n.constructs, ncol = n.constructs) - 
                       matrix(normalized.self, nrow = n.constructs, ncol = n.constructs)
  weight.matrix <- num.weight.matrix / den.weigth.matrix
  
  # Build wimp object
  wimp <- list()
  class(wimp) <- c("wimp", "list")

  # Global: metadata
  global_list$scale <- c(scale.min, scale.max)
  global_list$n.constructs <- as.numeric(n.constructs)
  global_list$wmatrix <- weight.matrix

  # Congruency classification
  congruency <- character(n.constructs)
  for(k in seq_len(n.constructs)){
    if(is.na(normalized.ideal[k]) || normalized.ideal[k] == 0){
      congruency[k] <- "dilemmatic"
    } else if(is.na(normalized.self[k]) || normalized.self[k] == 0){
      congruency[k] <- "undefined"
    } else if(sign(normalized.self[k]) == sign(normalized.ideal[k])){
      congruency[k] <- "congruent"
    } else {
      congruency[k] <- "discrepant"
    }
  }
  
  # Vertices: one row per construct
  vertices_df <- data.frame(
    id = seq_len(n.constructs),
    lpole = left.poles,
    rpole = right.poles,
    self = normalized.self,
    ideal = normalized.ideal,
    self_pole = self.poles,
    ideal_pole = ideal.poles,
    congruency = congruency,
    stringsAsFactors = FALSE
  )
  
  # Append extra attributes to vertices
  if(!is.null(extra_attrs) && ncol(extra_attrs) > 0){
    for(col_idx in seq_len(ncol(extra_attrs))){
      col_name <- if(!is.null(extra_attr_names) && col_idx <= length(extra_attr_names)) 
                    extra_attr_names[col_idx] else paste0("attr_", col_idx)
      if(is.na(col_name) || col_name == "") col_name <- paste0("attr_", col_idx)
      vertices_df[[col_name]] <- as.character(extra_attrs[[col_idx]])
    }
  }
  
  # Edges: one row per non-zero weight (excluding self-loops)
  edges_list <- vector("list", 0)
  for(i in seq_len(n.constructs)){
    for(j in seq_len(n.constructs)){
      if(i == j) next
      wval <- as.numeric(weight.matrix[i, j])
      if(is.na(wval) || abs(wval) < .Machine$double.eps) next
      edges_list[[length(edges_list) + 1]] <- list(
        id = paste(i, "t", j, sep = ""),
        from = as.integer(i),
        to = as.integer(j),
        weight = wval
      )
    }
  }
  
  edges_df <- if(length(edges_list) > 0){
    tmp <- do.call(rbind, lapply(edges_list, function(z) as.data.frame(z, stringsAsFactors = FALSE)))
    tmp$from <- as.integer(tmp$from)
    tmp$to <- as.integer(tmp$to)
    tmp$weight <- as.numeric(tmp$weight)
    tmp
  } else {
    data.frame(id = character(), from = integer(), to = integer(), weight = numeric(), stringsAsFactors = FALSE)
  }
  
  wimp$global <- global_list
  wimp$vertices <- vertices_df
  wimp$edges <- edges_df
  
  return(wimp)
}
