## RepGrid functions (Repertory Grid, Kelly 1955)
# Internal plotly implementations that keep the WimpTools look and feel.
# The biplot follows the approach of OpenRepGrid::biplot2d() (Heckmann, 2012):
# SVD of the centred ratings matrix with row/column scaling exponents g and h.

# Internal: extract ratings, construct poles and element names -------------------

.rg_extract <- function(x) {
  if (is.matrix(x) || is.data.frame(x)) {
    m <- as.matrix(x)
    if (!is.numeric(m)) stop("A ratings matrix must be numeric.")
    rn <- rownames(m); if (is.null(rn)) rn <- paste0("C", seq_len(nrow(m)))
    cn <- colnames(m); if (is.null(cn)) cn <- paste0("E", seq_len(ncol(m)))
    poles <- strsplit(rn, "\\s*(-|\\|)\\s*")
    left  <- vapply(poles, function(p) p[1], "")
    right <- vapply(poles, function(p) if (length(p) > 1) p[2] else "", "")
    return(list(ratings = m, left = left, right = right, elements = cn))
  }
  if (!methods::is(x, "repgrid")) {
    stop("`x` must be an OpenRepGrid 'repgrid' object or a numeric matrix.")
  }
  if (!requireNamespace("OpenRepGrid", quietly = TRUE)) {
    stop("Package 'OpenRepGrid' is required to use 'repgrid' objects.")
  }
  m  <- unname(OpenRepGrid::ratings(x))
  cn <- OpenRepGrid::constructs(x)
  list(
    ratings  = m,
    left     = as.character(cn$leftpole),
    right    = as.character(cn$rightpole),
    elements = OpenRepGrid::elements(x)
  )
}

# Internal: SVD-based biplot coordinates ---------------------------------------

.rg_biplot_coords <- function(rg, dims = 2, center = c("constructs", "none",
                              "elements", "both"), g = 0, h = 1 - g) {
  center <- match.arg(center)
  X <- rg$ratings
  if (anyNA(X)) stop("Ratings contain missing values.")
  X <- switch(center,
    constructs = X - rowMeans(X),
    elements   = sweep(X, 2, colMeans(X)),
    both       = { X <- X - rowMeans(X); sweep(X, 2, colMeans(X)) },
    none       = X
  )
  if (min(dim(X)) < dims) stop("Grid is too small for ", dims, " dimensions.")
  s <- svd(X)
  d <- s$d[seq_len(dims)]
  U <- s$u[, seq_len(dims), drop = FALSE]
  V <- s$v[, seq_len(dims), drop = FALSE]
  constructs <- sweep(U, 2, d^g, `*`)
  elements   <- sweep(V, 2, d^h, `*`)
  var_exp <- 100 * s$d^2 / sum(s$d^2)
  list(constructs = constructs, elements = elements,
       var_exp = var_exp[seq_len(dims)])
}

# Internal: 8-direction text anchor from a point's angle to the origin, so
# labels fan out away from their own vector instead of all stacking at a
# fixed "top center" position. Without this, biplots with many constructs
# become unreadable label soup as soon as the plot is narrower than ~900px
# (e.g. a dashboard widget card on a 1080p screen or smaller).
.rg_textpos <- function(x, y) {
  ang <- (atan2(y, x) * 180 / pi + 360) %% 360
  dirs <- c("middle right", "top right", "top center", "top left",
            "middle left", "bottom left", "bottom center", "bottom right")
  idx <- floor(((ang + 22.5) %% 360) / 45) + 1
  dirs[idx]
}

# repgrid_biplot ---------------------------------------------------------------

#' RepGrid Biplot - repgrid_biplot()
#'
#' @description Interactive plotly biplot (2D or 3D) of a repertory grid.
#'   Elements are plotted as labelled points and constructs as vectors
#'   spanning both poles.
#'
#' @param x An \code{OpenRepGrid} \code{repgrid} object, or a numeric matrix
#'   (constructs in rows, elements in columns; row names as
#'   \code{"left pole - right pole"}).
#' @param dim \code{2} (default) or \code{3}.
#' @param center Centring of the ratings before the SVD: \code{"constructs"}
#'   (default), \code{"none"}, \code{"elements"} or \code{"both"}.
#' @param g,h Exponents of the singular values applied to constructs and
#'   elements. Default \code{g = 0}, \code{h = 1 - g}.
#' @param text_size Text size multiplier. Defaults to 1.
#'
#' @return A plotly object.
#' @author Alejandro Sanfeliciano
#' @export
#' @examples
#' set.seed(1)
#' m <- matrix(sample(1:5, 40, replace = TRUE), 8,
#'             dimnames = list(paste0("left", 1:8, " - right", 1:8),
#'                             paste0("E", 1:5)))
#' repgrid_biplot(m)
#' repgrid_biplot(m, dim = 3)
repgrid_biplot <- function(x, dim = 2, center = "constructs", g = 0,
                           h = 1 - g, text_size = 1) {
  if (!dim %in% c(2, 3)) stop("`dim` must be 2 or 3.")
  rg <- .rg_extract(x)
  co <- .rg_biplot_coords(rg, dims = dim, center = center, g = g, h = h)

  # Scale construct vectors so they fill the element cloud
  el <- co$elements
  cs <- co$constructs
  k  <- max(abs(el)) / max(abs(cs)) * 0.9
  cs <- cs * k

  axis_lab <- function(i) {
    sprintf("<b>PC%d</b> [%.1f%%]", i, co$var_exp[i])
  }
  # Vector colours: pole of the ideal in green, the opposite pole in red and
  # dilemmatic constructs (ideal at the scale mid-point) in mustard
  col_ideal <- "#8cc63f"; col_opp <- "#F52722"; col_dil <- "#D4A017"
  col_none  <- "#9e9e9e"
  ps <- .rg_pole_status(x, rg)
  half_col <- function(i, side) {
    if (ps$dilemmatic[i]) return(col_dil)
    if (is.na(ps$ideal_side[i])) return(col_none)
    if (ps$ideal_side[i] == side) col_ideal else col_opp
  }
  col_r <- vapply(seq_len(nrow(cs)), half_col, "", side = "right")
  col_l <- vapply(seq_len(nrow(cs)), half_col, "", side = "left")
  # First column = actual self (blue), last column = ideal self (green)
  n_el    <- nrow(el)
  el_col  <- rep("#9e9e9e", n_el)
  el_col[1] <- "#2b6cb0"
  el_col[n_el] <- "#8cc63f"
  el_txt  <- ifelse(el_col == "#9e9e9e", "#555555", el_col)
  lab_l <- ifelse(nzchar(rg$left),  rg$left,  "")
  lab_r <- ifelse(nzchar(rg$right), rg$right, "")
  short <- function(v, n = 22) ifelse(nchar(v) > n, paste0(substr(v, 1, n - 1), "\u2026"), v)
  el_lab <- short(rg$elements)
  el_lab[c(1, n_el)] <- paste0("<b>", el_lab[c(1, n_el)], "</b>")  # self and ideal in bold
  hover_c <- paste0(lab_l, ifelse(nzchar(lab_r), " — ", ""), lab_r)

  if (dim == 2) {
    rng <- max(abs(rbind(el, cs, -cs))) * 1.15
    fig <- plotly::plot_ly(type = "scatter", mode = "markers")
    for (i in seq_len(nrow(cs))) {
      fig <- plotly::add_trace(
        fig, x = c(0, cs[i, 1]), y = c(0, cs[i, 2]),
        type = "scatter", mode = "lines", hoverinfo = "text",
        text = hover_c[i], showlegend = FALSE,
        line = list(color = col_r[i], dash = "dot", width = 1.25)
      )
      fig <- plotly::add_trace(
        fig, x = c(0, -cs[i, 1]), y = c(0, -cs[i, 2]),
        type = "scatter", mode = "lines", hoverinfo = "text",
        text = hover_c[i], showlegend = FALSE,
        line = list(color = col_l[i], dash = "dot", width = 1.25)
      )
    }
    # Fan each label out along its own vector/point angle instead of
    # stacking every one of them at a fixed "top center" position.
    pos_r <- .rg_textpos(cs[, 1], cs[, 2])
    pos_l <- .rg_textpos(-cs[, 1], -cs[, 2])
    pos_e <- .rg_textpos(el[, 1], el[, 2])
    fig <- fig %>%
      plotly::add_trace(
        x = cs[, 1], y = cs[, 2], type = "scatter", mode = "text", text = lab_r,
        textposition = pos_r, hoverinfo = "none", showlegend = FALSE,
        textfont = list(size = 12 * text_size, color = col_r)) %>%
      plotly::add_trace(
        x = -cs[, 1], y = -cs[, 2], type = "scatter", mode = "text", text = lab_l,
        textposition = pos_l, hoverinfo = "none", showlegend = FALSE,
        textfont = list(size = 12 * text_size, color = col_l)) %>%
      plotly::add_trace(
        x = el[, 1], y = el[, 2], type = "scatter", mode = "markers+text",
        text = el_lab, hovertext = rg$elements, textposition = pos_e, hoverinfo = "text",
        marker = list(color = el_col, size = 9),
        textfont = list(size = 15 * text_size, color = el_txt),
        showlegend = FALSE) %>%
      plotly::layout(
        xaxis = list(title = axis_lab(1), range = c(-rng, rng),
                     zeroline = TRUE, zerolinecolor = "black",
                     zerolinewidth = 2, showline = FALSE),
        yaxis = list(title = axis_lab(2), range = c(-rng, rng),
                     zeroline = TRUE, zerolinecolor = "black",
                     zerolinewidth = 2, showline = FALSE),
        margin = list(r = 70),
        showlegend = FALSE)
    return(fig)
  }

  # 3D ---------------------------------------------------------------------
  rng <- max(abs(rbind(el, cs))) * 1.15
  fig <- plotly::plot_ly()
  for (i in seq_len(nrow(cs))) {
    fig <- plotly::add_trace(
      fig, x = c(0, cs[i, 1]), y = c(0, cs[i, 2]), z = c(0, cs[i, 3]),
      type = "scatter3d", mode = "lines", hoverinfo = "text",
      text = hover_c[i], showlegend = FALSE,
      line = list(color = col_r[i], dash = "dot", width = 3)
    )
    fig <- plotly::add_trace(
      fig, x = c(0, -cs[i, 1]), y = c(0, -cs[i, 2]), z = c(0, -cs[i, 3]),
      type = "scatter3d", mode = "lines", hoverinfo = "text",
      text = hover_c[i], showlegend = FALSE,
      line = list(color = col_l[i], dash = "dot", width = 3)
    )
  }
  fig %>%
    plotly::add_trace(
      x = cs[, 1], y = cs[, 2], z = cs[, 3], type = "scatter3d",
      mode = "text", text = lab_r, hoverinfo = "none", showlegend = FALSE,
      textfont = list(size = 11 * text_size, color = col_r)) %>%
    plotly::add_trace(
      x = -cs[, 1], y = -cs[, 2], z = -cs[, 3], type = "scatter3d",
      mode = "text", text = lab_l, hoverinfo = "none", showlegend = FALSE,
      textfont = list(size = 11 * text_size, color = col_l)) %>%
    plotly::add_trace(
      x = el[, 1], y = el[, 2], z = el[, 3], type = "scatter3d",
      mode = "markers+text", text = el_lab, hovertext = rg$elements, hoverinfo = "text",
      textposition = "top center", showlegend = FALSE,
      marker = list(color = el_col, size = 5),
      textfont = list(size = 13 * text_size, color = el_txt)) %>%
    plotly::layout(scene = list(
      xaxis = list(title = axis_lab(1), range = c(-rng, rng)),
      yaxis = list(title = axis_lab(2), range = c(-rng, rng)),
      zaxis = list(title = axis_lab(3), range = c(-rng, rng)),
      aspectmode = "cube"), margin = list(r = 70))
}


# Internal: ideal pole and dilemmatic flag of each construct ---------------------
# Dilemmatic constructs are those where the ideal is rated at the scale mid-point.
.rg_pole_status <- function(x, rg) {
  M <- rg$ratings
  sc <- if (methods::is(x, "repgrid")) OpenRepGrid::getScale(x) else range(M, na.rm = TRUE)
  mid <- mean(as.numeric(sc))
  ideal <- M[, ncol(M)]
  side <- ifelse(ideal > mid, "right", ifelse(ideal < mid, "left", NA_character_))
  list(ideal_side = side, dilemmatic = ideal == mid)
}

# Internal: label shortener ----------------------------------------------------
.rg_short <- function(v, n = 35) {
  ifelse(nchar(v) > n, paste0(substr(v, 1, n - 1), "\u2026"), v)
}

# Internal: hclust + dendrogram segment layout (root left, leaves right) -------
.rg_hc_layout <- function(M, dist, method) {
  n <- nrow(M)
  hc <- stats::hclust(stats::dist(M, method = dist), method = method)
  pos <- numeric(n); pos[hc$order] <- seq_len(n)
  node_pos <- numeric(n - 1)
  segs <- vector("list", n - 1)
  get_pos <- function(k) if (k < 0) pos[-k] else node_pos[k]
  get_h   <- function(k) if (k < 0) 0 else hc$height[k]
  for (i in seq_len(n - 1)) {
    a <- hc$merge[i, 1]; b <- hc$merge[i, 2]
    pa <- get_pos(a); pb <- get_pos(b)
    ha <- get_h(a);   hb <- get_h(b);  h <- hc$height[i]
    node_pos[i] <- (pa + pb) / 2
    segs[[i]] <- list(x = c(ha, h, h, hb), y = c(pa, pa, pb, pb))
  }
  list(hc = hc, pos = pos, node_pos = node_pos,
       sx = unlist(lapply(segs, function(s) c(s$x, NA))),
       sy = unlist(lapply(segs, function(s) c(s$y, NA))))
}

# Internal: labels and leaf colours of the cluster dendrograms -----------------
# constructs: leaf colour by congruence (green congruent, red discrepant,
#   mustard dilemmatic = ideal at the mid-point, grey undefined) and the pole
#   associated with the self in bold.
# elements: actual self (first) in blue and ideal self (last) in green, both bold.
.rg_cluster_labels <- function(x, rg, along) {
  if (along == "elements") {
    n <- length(rg$elements)
    plain <- rg$elements
    tick <- .rg_short(plain)
    col <- rep("#555555", n)
    col[1] <- "#2b6cb0"; col[n] <- "#8cc63f"
    tick[c(1, n)] <- sprintf("<b><span style='color:%s'>%s</span></b>", col[c(1, n)], tick[c(1, n)])
    return(list(plain = plain, tick = tick, col = col))
  }
  M <- rg$ratings
  n <- nrow(M)
  sc <- if (methods::is(x, "repgrid")) OpenRepGrid::getScale(x) else range(M, na.rm = TRUE)
  mid <- mean(as.numeric(sc))
  self_side <- ifelse(M[, 1] < mid, "left", ifelse(M[, 1] > mid, "right", NA_character_))
  cls <- rep("neither", n)
  if (methods::is(x, "repgrid")) {
    cc <- tryCatch(OpenRepGrid::indexDilemma(x)$construct_classification$Classification,
                   error = function(e) NULL)
    if (length(cc) == n) cls <- as.character(cc)
  }
  cls[M[, ncol(M)] == mid] <- "dilemmatic"
  col <- c(congruent = "#8cc63f", discrepant = "#F52722",
           dilemmatic = "#D4A017", neither = "#9e9e9e")[cls]
  l <- .rg_short(rg$left, 20); r <- .rg_short(rg$right, 20)
  hasr <- nzchar(rg$right)
  lb <- ifelse(!is.na(self_side) & self_side == "left",  paste0("<b>", l, "</b>"), l)
  rb <- ifelse(!is.na(self_side) & self_side == "right", paste0("<b>", r, "</b>"), r)
  tick <- ifelse(hasr, paste0(lb, " \u2014 ", rb), lb)
  plain <- ifelse(hasr, paste0(rg$left, " \u2014 ", rg$right), rg$left)
  list(plain = plain, tick = tick, col = unname(col))
}

# Internal: precompute every dist x method layout for the widget menu ----------
.rg_cluster_options <- function(x, along, dists, methods) {
  rg <- .rg_extract(x)
  M <- if (along == "elements") t(rg$ratings) else rg$ratings
  lab <- .rg_cluster_labels(x, rg, along)$tick
  out <- list()
  for (d in dists) for (m in methods) {
    ly <- .rg_hc_layout(M, d, m)
    out[[paste(d, m, sep = "|")]] <- list(
      sx = ly$sx, sy = ly$sy, pos = ly$pos,
      hx = ly$hc$height, hy = ly$node_pos,
      ticktext = lab[ly$hc$order],
      hmax = max(ly$hc$height))
  }
  out
}

.rg_cluster_dists   <- c("euclidean", "manhattan", "maximum", "canberra")
.rg_cluster_methods <- c("ward.D", "ward.D2", "single", "complete", "average", "mcquitty")

# repgrid_cluster --------------------------------------------------------------

#' RepGrid Cluster Analysis - repgrid_cluster()
#'
#' @description Interactive plotly dendrogram from a hierarchical cluster
#'   analysis of the constructs or the elements of a repertory grid.
#'
#' @param x An \code{OpenRepGrid} \code{repgrid} object, or a numeric matrix
#'   (see \code{\link{repgrid_biplot}}).
#' @param along \code{"constructs"} (default) or \code{"elements"}.
#' @param dist Distance measure passed to \code{\link[stats]{dist}}.
#'   Defaults to \code{"euclidean"}.
#' @param method Agglomeration method passed to \code{\link[stats]{hclust}}.
#'   Defaults to \code{"ward.D"}.
#' @param text_size Text size multiplier. Defaults to 1.
#'
#' @return A plotly object.
#' @author Alejandro Sanfeliciano
#' @importFrom stats dist hclust
#' @export
#' @examples
#' set.seed(1)
#' m <- matrix(sample(1:5, 40, replace = TRUE), 8,
#'             dimnames = list(paste0("left", 1:8, " - right", 1:8),
#'                             paste0("E", 1:5)))
#' repgrid_cluster(m)
#' repgrid_cluster(m, along = "elements")
repgrid_cluster <- function(x, along = c("constructs", "elements"),
                            dist = "euclidean", method = "ward.D",
                            text_size = 1) {
  along <- match.arg(along)
  rg <- .rg_extract(x)
  M <- rg$ratings
  if (along == "elements") M <- t(M)
  n <- nrow(M)
  if (n < 3) stop("At least 3 ", along, " are needed for a cluster analysis.")

  lb <- .rg_cluster_labels(x, rg, along)
  lab <- lb$plain
  lab_col <- lb$col

  ly <- .rg_hc_layout(M, dist, method)
  pos <- ly$pos; node_pos <- ly$node_pos; hc <- ly$hc
  sx <- ly$sx; sy <- ly$sy
  short <- .rg_short
  ord <- hc$order
  hmax <- max(hc$height)

  fig <- plotly::plot_ly() %>%
    plotly::add_trace(
      x = sx, y = sy, type = "scatter", mode = "lines", hoverinfo = "none",
      line = list(color = "#9a9a9a", width = 1.5), showlegend = FALSE) %>%
    plotly::add_trace(
      x = rep(0, n), y = pos, type = "scatter", mode = "markers",
      text = lab, hoverinfo = "text", showlegend = FALSE,
      marker = list(color = lab_col, size = 8)) %>%
    plotly::add_trace(
      x = hc$height, y = node_pos, type = "scatter", mode = "markers",
      text = sprintf("d = %.2f", hc$height), hoverinfo = "text",
      showlegend = FALSE, marker = list(color = "#9a9a9a", size = 4)) %>%
    plotly::layout(
      xaxis = list(title = "<b>Distance</b>", autorange = "reversed",
                   range = c(hmax * 1.05, -hmax * 0.02), zeroline = FALSE,
                   showgrid = TRUE),
      yaxis = list(side = "right", tickmode = "array", tickvals = seq_len(n),
                   ticktext = lb$tick[ord], tickfont = list(size = 12 * text_size),
                   autorange = "reversed", showgrid = FALSE, zeroline = FALSE,
                   title = "", automargin = FALSE),
      margin = list(l = 40, r = 330),
      showlegend = FALSE)
  fig
}

# repgrid_dilemmas -------------------------------------------------------------

# Internal: implicative dilemma results via OpenRepGrid::indexDilemma() --------
.rg_dilemma_data <- function(x, ...) {
  if (!methods::is(x, "repgrid")) {
    stop("Implicative dilemmas need an OpenRepGrid 'repgrid' object.")
  }
  if (!requireNamespace("OpenRepGrid", quietly = TRUE)) {
    stop("Package 'OpenRepGrid' is required.")
  }
  d <- OpenRepGrid::indexDilemma(x, ...)
  cc <- d$construct_classification
  ga <- OpenRepGrid::constructs(d$grid_aligned)  # aligned: self pole on the left
  df <- d$dilemmas_df
  list(
    classification = data.frame(
      construct = as.character(cc$Construct),
      left = as.character(ga$leftpole), right = as.character(ga$rightpole),
      self = cc$Self, ideal = cc$Ideal,
      difference = cc$Difference, class = as.character(cc$Classification),
      stringsAsFactors = FALSE),
    dilemmas = data.frame(
      id_c = as.integer(df$id_c), congruent = sub("^[0-9]+\\. ", "", df$Congruent),
      id_d = as.integer(df$id_d), discrepant = sub("^[0-9]+\\. ", "", df$Discrepant),
      r = as.numeric(df$R), stringsAsFactors = FALSE),
    n_ids = d$no_ids, pid = d$measures$pid, iid = d$measures$iid,
    picid = d$measures$picid, r_min = d$r.min
  )
}

#' RepGrid Implicative Dilemmas Table - repgrid_dilemmas_table()
#'
#' @description Table of the implicative dilemmas (congruent and discrepant
#'   construct pairs correlated above a threshold) of a repertory grid.
#'
#' @param x An \code{OpenRepGrid} \code{repgrid} object. The self is assumed
#'   in the first column and the ideal in the last one.
#' @param ... Additional arguments passed to
#'   \code{OpenRepGrid::indexDilemma}.
#'
#' @return A data frame with one row per implicative dilemma.
#' @author Alejandro Sanfeliciano
#' @export
repgrid_dilemmas_table <- function(x, ...) {
  .rg_dilemma_data(x, ...)$dilemmas
}

#' RepGrid Implicative Dilemmas Plot - repgrid_dilemmas()
#'
#' @description Interactive plotly diagram of the implicative dilemmas:
#'   congruent constructs on the left, discrepant constructs on the right and
#'   a line for each dilemma whose width grows with the correlation.
#'   Constructs not involved in any dilemma are shown faded.
#'
#' @param x An \code{OpenRepGrid} \code{repgrid} object. The self is assumed
#'   in the first column and the ideal in the last one.
#' @param text_size Text size multiplier. Defaults to 1.
#' @param only_involved If \code{TRUE}, only the constructs that take part in
#'   at least one dilemma are drawn. Defaults to \code{FALSE}.
#' @param ... Additional arguments passed to
#'   \code{OpenRepGrid::indexDilemma}.
#'
#' @return A plotly object.
#' @author Alejandro Sanfeliciano
#' @export
repgrid_dilemmas <- function(x, text_size = 1, only_involved = FALSE, ...) {
  dd <- .rg_dilemma_data(x, ...)
  cl <- dd$classification
  dl <- dd$dilemmas
  cong <- which(cl$class == "congruent")
  disc <- which(cl$class == "discrepant")
  if (only_involved) {
    cong <- cong[cong %in% dl$id_c]
    disc <- disc[disc %in% dl$id_d]
  }
  if (length(cong) == 0 || length(disc) == 0) {
    msg <- if (only_involved) "No implicative dilemmas" else "No congruent / discrepant constructs"
    return(plotly::plot_ly(type = "scatter", mode = "markers") %>%
      plotly::layout(
        xaxis = list(visible = FALSE), yaxis = list(visible = FALSE),
        annotations = list(list(
          text = msg, showarrow = FALSE,
          font = list(size = 15, color = "#999999")))))
  }

  # Self pole (left, after alignment) in bold. Kept short (16 chars) because
  # the label sits in a fixed-width margin outside the plotting area (see
  # `margin` below) - a longer truncation needs a wider margin than most
  # dashboard widget cards have on a 1080p screen or smaller, which cuts the
  # text off against the left/right edge instead of just wrapping it.
  lab_html <- function(k) {
    l <- .rg_short(cl$left[k], 16); r <- .rg_short(cl$right[k], 16)
    paste0("<b>", l, "</b> - ", r)
  }
  # Both columns hang from the top; any blank space is left at the bottom
  top <- max(length(cong), length(disc))
  ypos <- function(ids) setNames(top + 1 - seq_along(ids), ids)
  yc <- ypos(cong); yd <- ypos(disc)
  col_c <- "#8cc63f"; col_d <- "#F52722"; col_e <- "#9a9a9a"

  inv_c <- cong %in% dl$id_c; inv_d <- disc %in% dl$id_d
  fig <- plotly::plot_ly(type = "scatter", mode = "markers")

  # Discrepant constructs on the left, congruent on the right
  for (i in seq_len(nrow(dl))) {
    y0 <- yd[as.character(dl$id_d[i])]; y1 <- yc[as.character(dl$id_c[i])]
    fig <- plotly::add_trace(
      fig, x = c(0, 1), y = c(y0, y1), type = "scatter", mode = "lines",
      hoverinfo = "none", showlegend = FALSE,
      line = list(color = col_e, width = 1 + 1.5 * dl$r[i]))
  }
  fig <- fig %>%
    plotly::add_trace(
      x = rep(0, length(disc)), y = unname(yd), type = "scatter",
      mode = "markers", hoverinfo = "text", hovertext = cl$construct[disc],
      showlegend = FALSE,
      marker = list(color = col_d, size = 11, opacity = ifelse(inv_d, 1, 0.35))) %>%
    plotly::add_trace(
      x = rep(1, length(cong)), y = unname(yc), type = "scatter",
      mode = "markers", hoverinfo = "text", hovertext = cl$construct[cong],
      showlegend = FALSE,
      marker = list(color = col_c, size = 11, opacity = ifelse(inv_c, 1, 0.35)))

  lab_ann <- function(x, y, text, inv, anchor, shift) list(
    x = x, y = y, xref = "x", yref = "y", text = text, showarrow = FALSE,
    xanchor = anchor, xshift = shift,
    font = list(size = 12 * text_size, color = if (inv) "#333333" else "#aaaaaa"))
  ann <- c(
    lapply(seq_along(disc), function(k)
      lab_ann(0, unname(yd)[k], lab_html(disc[k]), inv_d[k], "right", -25)),
    lapply(seq_along(cong), function(k)
      lab_ann(1, unname(yc)[k], lab_html(cong[k]), inv_c[k], "left", 25))
  )
  # Correlation values (.xx) on a white box that hides the line behind them
  val_ann <- lapply(seq_len(nrow(dl)), function(i) list(
    x = 0.5, y = (yd[as.character(dl$id_d[i])] + yc[as.character(dl$id_c[i])]) / 2,
    xref = "x", yref = "y", text = sub("^0", "", sprintf("%.2f", dl$r[i])),
    showarrow = FALSE, bgcolor = "#ffffff", borderpad = 2,
    font = list(size = 12 * text_size, color = "#666666")))
  ann <- c(ann, val_ann)

  fig %>% plotly::layout(
    xaxis = list(visible = FALSE, range = c(-0.05, 1.05), fixedrange = TRUE),
    yaxis = list(visible = FALSE, range = c(0.3, max(length(cong), length(disc)) + 0.7),
                 fixedrange = TRUE),
    # Was l = 310, r = 350: sized for full-width exports only. On a narrower
    # widget card (e.g. a dashboard grid cell on a 1080p screen or smaller)
    # that left almost no room for the actual plot area, cutting the pole
    # labels off against the edges. Shortened labels above let this shrink
    # to a still-generous but no longer overflowing margin.
    margin = list(l = 220, r = 220, t = 20, b = 20),
    annotations = ann, showlegend = FALSE)
}

# repgrid_indices --------------------------------------------------------------

#' RepGrid Cognitive Indices - repgrid_indices()
#'
#' @description Computes the main cognitive indices of a repertory grid with
#'   \code{OpenRepGrid}: intensity, PVAFF, Bieri, variability, polarization,
#'   bias, three conflict measures, self-construction distances and the
#'   dilemmatic constructs index. The self is assumed in the
#'   first column and the ideal in the last one.
#'
#' @param x An \code{OpenRepGrid} \code{repgrid} object.
#'
#' @return A list with three data frames: \code{global} (one row per index,
#'   with columns \code{group}, \code{key}, \code{value} and \code{unit}),
#'   \code{constructs} and \code{elements} (per construct / element intensity,
#'   polarization and conflict percentage).
#' @author Alejandro Sanfeliciano
#' @export
repgrid_indices <- function(x) {
  if (!methods::is(x, "repgrid")) stop("`x` must be an OpenRepGrid 'repgrid' object.")
  if (!requireNamespace("OpenRepGrid", quietly = TRUE)) stop("Package 'OpenRepGrid' is required.")
  ne <- ncol(OpenRepGrid::ratings(x))
  nc <- nrow(OpenRepGrid::ratings(x))
  try_ <- function(expr) tryCatch(suppressMessages(suppressWarnings(expr)), error = function(e) NULL)
  num <- function(v) if (is.null(v) || length(v) == 0) NA_real_ else as.numeric(v)[1]

  int  <- try_(OpenRepGrid::indexIntensity(x))
  pol  <- try_(OpenRepGrid::indexPolarization(x))
  c1   <- try_(OpenRepGrid::indexConflict1(x))
  c2   <- try_(OpenRepGrid::indexConflict2(x))
  c3   <- try_(OpenRepGrid::indexConflict3(x))
  sc   <- try_(OpenRepGrid::indexSelfConstruction(x, self = 1, ideal = ne))
  dlm  <- try_(OpenRepGrid::indexDilemmatic(x, ideal = ne))

  row <- function(group, key, value, unit = "") data.frame(
    group = group, key = key, value = value, unit = unit, stringsAsFactors = FALSE)
  global <- rbind(
    row("structure", "int_total", num(int$total.int)),
    row("structure", "int_c", num(int$c.int.mean)),
    row("structure", "int_e", num(int$e.int.mean)),
    row("structure", "pvaff", 100 * num(try_(OpenRepGrid::indexPvaff(x))), "%"),
    row("structure", "bieri", num(try_(OpenRepGrid::indexBieri(x))$bieri)),
    row("structure", "variability", num(try_(OpenRepGrid::indexVariability(x)))),
    row("structure", "bias", num(try_(OpenRepGrid::indexBias(x)))),
    row("structure", "polarization", 100 * num(pol$polarization_total$Polarization), "%"),
    row("conflict", "conf1", 100 * num(c1$prop.imbalanced), "%"),
    row("conflict", "conf2", 100 * num(c2$prop.imbalanced), "%"),
    row("conflict", "conf3", num(c3$overall), "%"),
    row("self", "self_ideal", num(sc$self_ideal)),
    row("self", "self_others", num(sc$self_others)),
    row("self", "ideal_others", num(sc$ideal_others)),
    row("self", "dilemmatic", 100 * num(dlm$perc_dilemmatic), "%")
  )

  rg <- .rg_extract(x)
  cn <- ifelse(nzchar(rg$right), paste0(rg$left, " — ", rg$right), rg$left)
  constructs <- data.frame(
    construct = cn,
    intensity = if (!is.null(int)) as.numeric(int$c.int) else NA_real_,
    polarization = if (!is.null(pol)) 100 * pol$polarization_constructs$Polarization else NA_real_,
    conflict = if (!is.null(c3)) as.numeric(c3$c.perc[[1]]) else NA_real_,
    stringsAsFactors = FALSE)
  elements <- data.frame(
    element = rg$elements,
    intensity = if (!is.null(int)) as.numeric(int$e.int) else NA_real_,
    polarization = if (!is.null(pol)) 100 * pol$polarization_elements$Polarization else NA_real_,
    conflict = if (!is.null(c3)) as.numeric(c3$e.perc[[1]]) else NA_real_,
    stringsAsFactors = FALSE)
  list(global = global, constructs = constructs, elements = elements)
}

# importrepgrid_json ----------------------------------------------------------------

#' Import a Repertory Grid from JSON -- importrepgrid_json()
#'
#' @description Reads a repertory grid (RepGrid) stored in the JSON format
#'   used by the PsychLab web application and returns an
#'   \code{OpenRepGrid} \code{repgrid} object, ready to be used with
#'   \code{\link{repgrid_biplot}}, \code{\link{repgrid_cluster}},
#'   \code{\link{repgrid_dilemmas}}, \code{\link{repgrid_indices}} and the
#'   \code{widget_repgrid_*} functions.
#'
#' @param x A path to a \code{.json} file, a JSON string, or a list already
#'   parsed with \code{jsonlite::fromJSON(simplifyVector = FALSE)}.
#'
#' @details The record must have \code{type = "repgrid"} (when present) and a
#'   \code{data} object with \code{scaleMin}, \code{scaleMax},
#'   \code{elements} (names) and \code{constructs}, each one with
#'   \code{left}, \code{right} and \code{ratings} (one value per element,
#'   \code{null} for a missing rating). By convention of the package the self
#'   is the first element and the ideal self the last one.
#'
#'   The rest of the record (\code{id}, \code{title}, \code{status},
#'   \code{patientId}, \code{patientName}, \code{notes}, \code{params},
#'   \code{createdAt}, \code{updatedAt}) is kept in the \code{meta} slot of the
#'   returned object, e.g. \code{x@meta$title}.
#'
#' @return An \code{OpenRepGrid} \code{repgrid} object.
#' @author Alejandro Sanfeliciano
#' @export
#' @examples
#' json <- '{"type": "repgrid", "title": "Demo", "data": {
#'   "scaleMin": 1, "scaleMax": 5,
#'   "elements": ["Self", "Mother", "Ideal"],
#'   "constructs": [
#'     {"left": "calm", "right": "anxious", "ratings": [2, 4, 1]},
#'     {"left": "open", "right": "closed", "ratings": [3, 5, 2]},
#'     {"left": "active", "right": "passive", "ratings": [1, 2, 1]}]}}'
#' rg <- importrepgrid_json(json)
#' rg@meta$title
importrepgrid_json <- function(x) {
  if (!requireNamespace("OpenRepGrid", quietly = TRUE)) {
    stop("Package 'OpenRepGrid' is required to import a RepGrid.")
  }
  rec <- if (is.list(x)) {
    x
  } else if (is.character(x) && length(x) == 1) {
    tryCatch(jsonlite::fromJSON(x, simplifyVector = FALSE),
             error = function(e) stop("`x` is not a valid JSON file or string: ",
                                      conditionMessage(e), call. = FALSE))
  } else {
    stop("`x` must be a path, a JSON string or a parsed list.")
  }

  type <- rec$type
  if (!is.null(type) && !identical(tolower(type), "repgrid")) {
    stop("The record has type '", type, "', not 'repgrid'.",
         if (identical(tolower(type), "wimpgrid")) " Use importwimp() for WimpGrids." else "",
         call. = FALSE)
  }
  d <- rec$data
  if (is.null(d)) stop("The record has no `data` field.", call. = FALSE)

  # Elements: plain strings or objects with a name
  els <- vapply(d$elements, function(e) {
    if (is.list(e)) as.character(e$name %||% e$label %||% "") else as.character(e)
  }, "")
  if (length(els) < 2) stop("At least 2 elements are needed.", call. = FALSE)
  cons <- d$constructs
  if (length(cons) < 2) stop("At least 2 constructs are needed.", call. = FALSE)

  ne <- length(els)
  rating_mat <- t(vapply(seq_along(cons), function(i) {
    r <- cons[[i]]$ratings
    if (length(r) != ne) {
      stop("Construct ", i, " has ", length(r), " ratings but there are ", ne,
           " elements.", call. = FALSE)
    }
    vapply(r, function(v) if (is.null(v)) NA_real_ else as.numeric(v), 0)
  }, numeric(ne)))
  if (anyNA(rating_mat)) {
    warning("The grid has missing ratings; some functions need complete grids.",
            call. = FALSE)
  }

  smin <- as.numeric(d$scaleMin %||% min(rating_mat, na.rm = TRUE))
  smax <- as.numeric(d$scaleMax %||% max(rating_mat, na.rm = TRUE))
  if (!is.finite(smin) || !is.finite(smax) || smin >= smax) {
    stop("`scaleMin` must be lower than `scaleMax`.", call. = FALSE)
  }
  if (any(rating_mat < smin | rating_mat > smax, na.rm = TRUE)) {
    stop("Some ratings are outside the scale [", smin, ", ", smax, "].", call. = FALSE)
  }

  pole <- function(k) vapply(cons, function(cn) as.character(cn[[k]] %||% ""), "")
  rg <- OpenRepGrid::makeRepgrid(list(
    name = els, l.name = pole("left"), r.name = pole("right"),
    scores = as.vector(t(rating_mat))))
  rg <- OpenRepGrid::setScale(rg, smin, smax)

  meta <- rec[setdiff(names(rec), "data")]
  rg@meta <- meta[!vapply(meta, is.null, TRUE)]
  rg
}

`%||%` <- function(a, b) if (is.null(a)) b else a
