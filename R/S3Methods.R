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
#' # print(example.wimp)
#'
#' @export
print.wimp <- function(x, n = 10, m = 10, digits = 3, ...){
  stopifnot(is.list(x), inherits(x, "wimp"))

  # Helpers ---------------------------------------------------------------
  .trunc <- function(z, width){
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
