## S3 methods for pretty console output

#' Print method for wimp objects
#'
#' Presents a concise, readable summary of a `wimp` object, including
#' global metadata, a preview of vertices, and a preview of edges.
#'
#' @param x A `wimp` object created by `importwimp()`.
#' @param n Number of vertex rows to show (head). Default: 10.
#' @param m Number of edge rows to show (head). Default: 10.
#' @param digits Number of digits for numeric rounding. Default: 3.
#' @param ... Ignored.
#'
#' @examples
#' # print(example_wimp)
#'
#' @export
print.wimp <- function(x, n = 10, m = 10, digits = 3, ...) {
  stopifnot(is.list(x), inherits(x, "wimp"))

  # Helpers ---------------------------------------------------------------
  .trunc <- function(z, width) {
    z <- as.character(z)
    too_long <- nchar(z) > width
    z[too_long] <- paste0(substr(z[too_long], 1, width - 1), "…")
    z
  }
  .num <- function(v) ifelse(is.na(v), NA, round(v, digits))
  .cat <- function(...) cat(paste0(...), sep = "")
  .rule <- function(title = NULL, width = getOption("width", 80)){
    if(is.null(title) || title == ""){
      .cat(paste0(strrep("-", width), "\n"))
    } else {
      ttl <- paste0(" ", title, " ")
      k <- max(0, width - nchar(ttl))
      left <- floor(k/2); right <- k - left
      .cat(strrep("-", left), ttl, strrep("-", right), "\n")
    }
  }
  .df_preview <- function(df, max_rows = 10, max_col_width = 24){
    if(NROW(df) == 0){ .cat("(empty)\n"); return(invisible()) }
    rows <- head(df, max_rows)

    # Determine alignment per column before coercion
    is_num_col <- vapply(rows, is.numeric, logical(1))

    # Format numbers, then coerce all to character and truncate
    for (nm in names(rows)){
      if(is.numeric(rows[[nm]])) rows[[nm]] <- .num(rows[[nm]])
    }
    rows[] <- lapply(rows, as.character)
    rows[] <- lapply(rows, .trunc, width = max_col_width)

    # Compute widths from header and data (after truncation)
    colnames_now <- names(rows)
    data_char <- as.data.frame(lapply(rows, as.character), stringsAsFactors = FALSE)
    width_data <- vapply(data_char, function(v){ if(length(v)==0) 0L else max(nchar(v), na.rm = TRUE) }, integer(1))
    width_head <- nchar(colnames_now)
    w <- pmax(width_data, width_head)
    w <- pmin(w, max_col_width)

    # Header
    hdr_cells <- mapply(function(h, ww){ format(h, width = ww, justify = "left") }, colnames_now, w, USE.NAMES = FALSE)
    hdr <- paste(hdr_cells, collapse = "  |  ")
    .cat(hdr, "\n")
    .cat(strrep("-", nchar(hdr)), "\n")

    # Rows with alignment: numeric right, others left
    mat <- as.matrix(data_char)
    for (i in seq_len(nrow(mat))){
      cells <- character(ncol(mat))
      for (j in seq_len(ncol(mat))){
        val <- mat[i, j]
        if(is.na(val)) val <- "NA"
        just <- if(is_num_col[j]) "right" else "left"
        cells[j] <- format(val, width = w[j], justify = just)
      }
      .cat(paste(cells, collapse = "  |  "), "\n")
    }
    if(NROW(df) > max_rows) .cat("… with ", NROW(df) - max_rows, " more rows\n")
  }

  width <- getOption("width", 80)

  # Header ----------------------------------------------------------------
  .rule("WimpGrid")

  # Global ----------------------------------------------------------------
  gl <- x$global
  n_con <- if(!is.null(gl$n.constructs)) as.integer(gl$n.constructs) else if(!is.null(x$vertices)) nrow(x$vertices) else NA_integer_
  scale <- if(!is.null(gl$scale)) paste0("[", paste(.num(gl$scale), collapse = ", "), "]") else "(not set)"

  # Weight matrix summary
  wdim <- "(none)"; edges_n <- NA_integer_; density <- NA_real_
  if(!is.null(gl$wmatrix) && is.matrix(gl$wmatrix)){
    d <- dim(gl$wmatrix)
    wdim <- paste0(d[1], " x ", d[2])
    nz <- sum(gl$wmatrix != 0, na.rm = TRUE) - sum(diag(gl$wmatrix) != 0, na.rm = TRUE)
    edges_n <- max(0L, nz)
    if(d[1] > 1) density <- edges_n / (d[1] * (d[2] - 1))
  }

  .cat("- Constructs: ", n_con, "\n")
  .cat("- Scale:     ", scale, "\n")
  .cat("- Matrix:    ", wdim)
  if(!is.na(edges_n)) .cat("  |  edges: ", edges_n)
  if(!is.na(density)) .cat("  |  density: ", sprintf("%.1f%%", 100 * density))
  .cat("\n")
  # Additional metadata (up to 6 key-value pairs, excluding known keys)
  if(is.list(gl)){
    meta_names <- setdiff(names(gl), c("scale","n.constructs","wmatrix"))
    if(length(meta_names) > 0){
      k <- head(meta_names, 6)
      .cat("- Metadata: ")
      kv <- vapply(k, function(nm){
        val <- gl[[nm]]
        if(length(val) > 1) val <- paste(val, collapse = ", ")
        paste0(nm, ": ", .trunc(val, 30))
      }, character(1))
      .cat(paste(kv, collapse = "  |  "), "\n")
      if(length(meta_names) > length(k)) .cat("            … and ", length(meta_names) - length(k), " more\n")
    }
  }
  .cat("\n")
  .cat("\n")

  # Vertices --------------------------------------------------------------
  .rule("Vertices", width)
  v <- x$vertices
  if(is.null(v) || !is.data.frame(v) || nrow(v) == 0){
    .cat("(no vertices)\n")
  } else {
    show_cols <- c("id","lpole","rpole","self","ideal","self_pole","ideal_pole","congruency")
    show_cols <- intersect(show_cols, names(v))
    # include common extra columns if present (first up to 3)
    extra <- setdiff(names(v), c(show_cols))
    pick_extra <- head(extra, 3)
    .df_preview(v[, c(show_cols, pick_extra), drop = FALSE], max_rows = n)
  }
  .cat("\n")
  .cat("\n")

  # Edges -----------------------------------------------------------------
  .rule("Edges", width)
  e <- x$edges
  if(is.null(e) || !is.data.frame(e) || nrow(e) == 0){
    # attempt to summarize from matrix
    if(!is.null(gl$wmatrix) && is.matrix(gl$wmatrix)){
      d <- dim(gl$wmatrix)
      nz <- sum(gl$wmatrix != 0, na.rm = TRUE) - sum(diag(gl$wmatrix) != 0, na.rm = TRUE)
      .cat("(no edges table) — inferred ", max(0L, nz), " non-zero weights from ", d[1], "x", d[2], " matrix\n")
    } else {
      .cat("(no edges)\n")
    }
  } else {
    # quick weight summary
    w <- e$weight
    if(is.numeric(w)){
      qs <- stats::quantile(w, probs = c(0, .25, .5, .75, 1), na.rm = TRUE)
      .cat("- Count: ", nrow(e), "  |  weight min/25%/50%/75%/max: ",
           paste(.num(qs), collapse = " / "), "\n")
    } else {
      .cat("- Count: ", nrow(e), "\n")
    }
    show_cols_e <- intersect(c("from","to","weight"), names(e))
    .df_preview(e[, show_cols_e, drop = FALSE], max_rows = m)
  }

  invisible(x)
}

#' Print method for scn objects
#'
#' Shows the parameters used to generate the scenario matrix and the values
#' for each iteration.
#'
#' @param x An object of class `scn` returned by `scenariomatrix()`.
#' @param digits Digits to round numeric values. Default: 3.
#' @param ... Ignored.
#'
#' @export
print.scn <- function(x, digits = 3, ...) {
  stopifnot(is.list(x), inherits(x, "scn"))

  `%||%` <- function(a, b) if (!is.null(a)) a else b

  params <- if (!is.null(x$params)) x$params else list()
  method <- if (!is.null(x$method)) x$method else list()

  act_vec <- params$act_vector
  act_vec_txt <- if (!is.null(act_vec)) {
    paste(round(act_vec, digits), collapse = ", ")
  } else {
    "(not provided)"
  }

  width <- getOption("width", 80)
  rule <- function(title = NULL) {
    if (is.null(title) || title == "") {
      cat(strrep("-", width), "\n", sep = "")
    } else {
      ttl <- paste0(" ", title, " ")
      k <- max(0, width - nchar(ttl))
      left <- floor(k / 2)
      right <- k - left
      cat(strrep("-", left), ttl, strrep("-", right), "\n", sep = "")
    }
  }

  rule("Scenario (scn)")
  cat("Parameters\n")
  cat("  inference:  ", `%||%`(params$infer, method$infer), "\n", sep = "")
  cat("  threshold:  ", `%||%`(params$threshold, method$threshold), "\n",
      sep = "")
  if (!is.null(params$max_iter)) {
    cat("  max_iter:   ", params$max_iter, "\n", sep = "")
  }
  if (!is.null(params$e)) {
    cat("  tolerance:  ", params$e, "\n", sep = "")
  }
  if (!is.null(params$stop_iter)) {
    cat("  stop_iter:  ", params$stop_iter, "\n", sep = "")
  }
  if (!is.null(params$exclude_dilemmatics)) {
    cat("  exclude_dilemmatics: ", params$exclude_dilemmatics, "\n",
        sep = "")
  }
  cat("  act_vector: ", act_vec_txt, "\n", sep = "")
  cat("  constructs: ", length(x$constructs$constructs), "\n", sep = "")
  cat("  iterations: ", nrow(x$values), "\n", sep = "")

  mat <- x$values
  rule("Values by iteration")
  if (is.null(mat) || length(mat) == 0) {
    cat("(no values)\n")
    return(invisible(x))
  }
  if (!is.matrix(mat)) {
    mat <- as.matrix(mat)
  }
  mat <- round(mat, digits)

  print(mat)

  invisible(x)
}

#' Print method for self_index results
#'
#' Provides a formatted display of self-ideal similarity index results including
#' global indices, construct classification, and detailed construct-level analysis.
#'
#' @param x A self_index object returned by the \code{\link{self_index}} function
#' @param ... Additional arguments (not used)
#' @return Invisibly returns the input object
#' @method print self_index
#' @export
print.self_index <- function(x, ...) {

  cat("\n")
  cat("===================================================================\n")
  cat("                     SELF ANALYSIS                              \n")
  cat("===================================================================\n")

  # Method information
  method_name <- switch(x$method[1],
                       "ssi" = "Similarity Self-Ideal  (SSI)",
                       "pearson" = "Pearson Correlation",
                       "spearman" = "Spearman Correlation",
                       "kendall" = "Kendall Correlation",
                       paste("Correlation -", x$method[1]))

  # RC correction status - SSI doesn't need RC correction
  if(x$method[1] == "ssi") {
    rc_status <- "Not applicable (SSI does not require RC correction)"
  } else {
    rc_status <- if(x$method[2] == "rc") "with RC Correction" else "without RC Correction"
  }

  cat("Method:     ", method_name, "\n")
  cat("Correction: ", rc_status, "\n")
  cat("\n")

  # Global indices
  cat("-------------------------------------------------------------------\n")
  cat("                        GLOBAL INDICES                        \n")
  cat("-------------------------------------------------------------------\n")

  global_data <- x$global
  cat(sprintf("%-20s %8.4f\n", "Self/Ideal:", global_data[["Self/Ideal"]]))
  cat(sprintf("%-20s %8.4f\n", "Self/Hypothetical:", global_data[["Self/Hypo"]]))
  cat(sprintf("%-20s %8.4f\n", "Ideal/Hypothetical:", global_data[["Ideal/Hypo"]]))

  cat("\n")

  # Construct level analysis
  cat("-------------------------------------------------------------------\n")
  cat("                    CONSTRUCT LEVEL ANALYSIS                  \n")
  cat("-------------------------------------------------------------------\n")

  construct_data <- x$construct

  # Get proper construct classification using construct_index function
  if(!is.null(x$wimp)) {
    ci_result <- construct_index(x$wimp)
    cat("\nConstruct Classification:\n")
    for(i in seq_len(nrow(ci_result))) {
      type <- rownames(ci_result)[i]
      count <- ci_result[i, "Frequency"]
      cat(sprintf("  %-12s: %2d constructs\n", type, count))
    }
  } else {
    # Fallback if wimp object not available
    congruence_summary <- table(construct_data[["Congruence Scenario"]])
    cat("\nConstruct Classification:\n")
    for(type in names(congruence_summary)) {
      cat(sprintf("  %-12s: %2d constructs\n", type, congruence_summary[type]))
    }
  }
  cat("\n")

  # Detailed table header
  cat(sprintf("%-35s %-12s %8s %8s\n",
              "Hypothetical Scenario",
              "Congruence",
              "Hypothetical/Self",
              "Hypothetical/Ideal"))
  cat("-------------------------------------------------------------------\n")

  # Print each construct
  for(i in seq_len(nrow(construct_data))) {
    scenario <- construct_data[[i, "Hypothetical Scenario"]]
    congruence <- construct_data[[i, "Congruence Scenario"]]
    self_sim <- construct_data[[i, "SHS"]]
    ideal_sim <- construct_data[[i, "SHI"]]

    # Truncate scenario name if too long
    if(nchar(scenario) > 32) {
      scenario <- paste0(substr(scenario, 1, 29), "...")
    }

    cat(sprintf("%-35s %-12s %8.4f %8.4f\n",
                scenario,
                congruence,
                self_sim,
                ideal_sim))
  }

  cat("\n")
  cat("===================================================================\n")
  cat("\n")

  invisible(x)
}

