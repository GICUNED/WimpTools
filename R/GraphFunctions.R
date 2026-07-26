## DIGRAPH FUNCTIONS ##

# Self Digraph -----------------------------------------------------------------

#' Self Digraph -- digraph()
#'
#' @description
#' Creates a directed graph (digraph) visualization representing the self of
#' the person being assessed based on their personal constructs and the
#' relationships between them. The digraph displays constructs as nodes with
#' colors indicating congruency between self and ideal values, and edges
#' representing the influence relationships derived from the weight matrix.
#'
#' @param wimp A subject's WimpGrid object. Must be a "wimp" S3 object
#'        imported by the \code{\link{importwimp}} function, or a "scn"
#'        scenario object from \code{\link{scenariomatrix}}.
#' @param vertex_vector Numeric vector defining the value of each vertex in the
#'        digraph. If \code{NA} (default), uses the normalized self values
#'        from the wimp object.
#' @param ideal_vector Numeric vector defining the ideal value of each vertex.
#'        If \code{NA} (default), uses the normalized ideal values from the
#'        wimp object.
#' @param width Character string specifying the graph width. Default is "100\%".
#' @param height Height of the digraph. Default is "90vh".
#' @param color Character string specifying the color palette. Options are
#'        "red/green" (default) and "grey scale".
#' @param layout Character string specifying the layout algorithm. Options are:
#'        "graphopt" (default), "circle", "rtcircle", "tree", "mds", "grid",
#'        and "areas".
#' @param show Logical vector or single value defining which constructs to
#'        display. Default is \code{TRUE} (show all constructs).
#' @param hide_direct Logical; if \code{TRUE}, hides positive direct
#'        relationships between nodes. Default is \code{FALSE}.
#' @param areas Logical; if \code{TRUE}, draws colored areas grouping nodes
#'        by the specified attribute. Default is \code{FALSE}.
#' @param area_attr Character string specifying the vertex attribute name for
#'        grouping nodes into areas. Default is "category". Must be a column
#'        name in the wimp vertices data frame.
#' @param area_color Character vector of hex colors for area backgrounds. If
#'        \code{NA} (default), uses predefined color palette.
#' @param pad_side Numeric value specifying padding in pixels for areas drawn
#'        around nodes. Default is 50.
#' @param rounding Numeric value specifying radius for rounding area corners
#'        (0 for sharp corners). Default is 10.
#' @param min_weight Numeric value specifying the minimum absolute weight for 
#'        an edge to be displayed. Default is 0 (show all edges).
#' @param interactive_options Logical; if \code{TRUE}, adds a comprehensive 
#'        options panel to the visualization for real-time adjustments of 
#'        palettes, layouts, and filters. Default is \code{TRUE}.
#'
#' @details
#' The digraph visualization provides insights into:
#' \itemize{
#'   \item \strong{Node colors}: Reflect congruency between self and ideal
#'   values
#'     \itemize{
#'       \item Red/green palette: Green (congruent), red (discrepant),
#'             grey (undefined), yellow (dilemmatic)
#'     }
#'   \item \strong{Node sizes}: Proportional to absolute self values
#'   \item \strong{Edge colors}: Represent relationship types
#'   (positive/negative)
#'   \item \strong{Edge directions}: Show influence flow between constructs
#' }
#'
#' @note
#' The function automatically reorients the weight matrix based on vertex signs
#' to ensure proper directionality. For reproducible layouts, the function uses
#' a fixed random seed (33) for graph layout algorithms.
#'
#' @author Alejandro Sanfeliciano
#'
#' @return A \code{visNetwork} interactive graph object that can be displayed
#'         in R Viewer, R Markdown documents, or Shiny applications.
#'
#' @seealso
#' \code{\link{idealdigraph}} for ideal self visualization,
#' \code{\link{simdigraph}} for scenario-based digraphs,
#' \code{\link{inout_digraph}} for construct relationship analysis
#'
#' @import visNetwork
#' @import plotly
#' @importFrom htmlwidgets JS onRender
#' @importFrom jsonlite toJSON
#' @importFrom magrittr %>%
#' @importFrom visNetwork visNetwork visOptions visInteraction visPhysics
#'             visIgraphLayout visEvents
#' @export
#'
#' @examples
#' # Basic digraph with default settings
#' digraph(example_wimp)
#'
#' # Custom layout and colors
#' digraph(example_wimp, layout = "circle", color = "grey scale")
#'
#' # Show areas grouped by category
#' digraph(example_wimp, areas = TRUE, area.attr = "category")

#'

digraph <- function(wimp, vertex_vector = NA, ideal_vector = NA, width = "100%",
                    height = "90vh", color = "red/green", layout = "graphopt",
                    show = TRUE, hide_direct = FALSE,
                    areas = FALSE, area_attr = "category", area_color = NA,
                    pad_side = 50, rounding = 10, min_weight = 0,
                    interactive_options = TRUE, sim_data = NULL, export_name = NULL, lang = "en", ...) {

  # --- Localization (see R/i18n.R) ---
  t <- wt_i18n(lang)

  if (is.null(export_name)) {
    export_name <- deparse(substitute(wimp))
  }

  # ==========================================
  # INPUT VALIDATION
  # ==========================================
  allowed_colors  <- c("red/green", "grey scale", "colorblind",
                       "pastel", "dark", "viridis")
  allowed_layouts <- c("graphopt", "circle", "rtcircle", "tree", "mds", "grid",
                       "areas")

  if (!(inherits(wimp, "wimp") || inherits(wimp, "scn"))) {
    stop("'wimp' must be an object of class 'wimp' or 'scn'.")
  }

  n_constructs <- if (inherits(wimp, "wimp")) {
    if (is.null(wimp$vertices) || !is.data.frame(wimp$vertices)) {
      stop("'wimp$vertices' must exist and be a data.frame in 'wimp' objects.")
    }
    nrow(wimp$vertices)
  } else {
    length(wimp$constructs[[1]])
  }

  if (!is.na(vertex_vector[1])) {
    if (!is.numeric(vertex_vector)) stop("'vertex_vector' must be numeric.")
    if (length(vertex_vector) != n_constructs) {
      stop("Length of 'vertex_vector' (", length(vertex_vector), ") does not 
           match number of constructs (", n_constructs, ").")
    }
  }

  if (!is.na(ideal_vector[1])) {
    if (!is.numeric(ideal_vector)) stop("'ideal_vector' must be numeric.")
    if (length(ideal_vector) != n_constructs) {
      stop("Length of 'ideal_vector' (", length(ideal_vector), ") does not match
           number of constructs (", n_constructs, ").")
    }
  }

  if (!(color %in% allowed_colors)) {
    stop("'color' must be one of: ", paste(allowed_colors, collapse = ", "))
  }
  if (!(layout %in% allowed_layouts)) {
    stop("'layout' must be one of: ", paste(allowed_layouts, collapse = ", "))
  }

  if (!is.logical(show)) {
    stop("'show' must be logical.")
  } else if (length(show) > 1 && length(show) != n_constructs) {
    stop("If 'show' is a vector it must have length equal to number of 
    constructs (", n_constructs, ").")
  }

  if (!is.logical(hide_direct) || length(hide_direct) != 1) stop("'hide_direct'
   must be logical of length 1.")
  if (!is.logical(areas) || length(areas) != 1) stop("'areas' must be logical of
   length 1.")

  if (!is.character(area_attr) || length(area_attr) != 1) stop("'area_attr' must
   be character of length 1.")
  if (areas && inherits(wimp, "wimp") &&
        !(area_attr %in% names(wimp$vertices))) {
    stop("'area_attr' (", area_attr, ") does not exist in 'wimp$vertices' 
         columns.")
  }

  if (!is.character(width) || length(width) != 1)  stop("'width' must be 
  character of length 1.")
  if (!is.character(height) || length(height) != 1)
    stop("'height' must be character of length 1.")

  if (!is.numeric(pad_side) || length(pad_side) != 1 || pad_side < 0) {
    stop("'pad_side' must be numeric >= 0.")
  }
  if (!is.numeric(rounding) || length(rounding) != 1 || rounding < 0) {
    stop("'rounding' must be numeric >= 0.")
  }
  if (!is.numeric(min_weight) || length(min_weight) != 1 || min_weight < 0) {
    stop("'min_weight' must be numeric >= 0.")
  }
  if (!is.logical(interactive_options) || length(interactive_options) != 1) {
    stop("'interactive_options' must be logical.")
  }

  if (inherits(wimp, "wimp") && !is.null(wimp$global$weight_matrix)) {
    wm <- wimp$global$weight_matrix
    if (!is.matrix(wm) || nrow(wm) != ncol(wm) || nrow(wm) != n_constructs) {
      stop("'wimp$global$weight_matrix' must be a square matrix with dimensions 
           equal to number of constructs.")
    }
  }

  # ==========================================
  # HELPER FUNCTIONS
  # ==========================================

  # Place nodes within a grid cell for areas layout
  .place_in_cell <- function(ids, cx, cy, cell_size) {
    n <- length(ids)
    if (n == 1L) {
      data.frame(id = ids, x = cx, y = cy)
    } else if (n == 2L) {
      r <- cell_size * 0.22
      ang <- c(pi / 6, 5 * pi / 6)
      data.frame(id = ids,
                 x = cx + r * cos(ang),
                 y = cy + r * sin(ang))
    } else {
      r <- cell_size * 0.28
      ang <- seq(0, 2 * pi, length.out = n + 1L)[- (n + 1L)]
      data.frame(id = ids,
                 x = cx + r * cos(ang),
                 y = cy + r * sin(ang))
    }
  }

  # Extract construct data from wimp object
  .extract_wimp_data <- function(wimp, area_attr) {
    if (inherits(wimp, "wimp")) {
      stopifnot(!is.null(wimp$vertices), is.data.frame(wimp$vertices))

      # Get weight matrix or create empty matrix
      wmatrix <- if (!is.null(wimp$global$weight_matrix)) {
        wimp$global$weight_matrix
      } else {
        matrix(0, nrow(wimp$vertices), nrow(wimp$vertices))
      }

      # Get area vector based on attribute
      area_vec <- if (area_attr %in% names(wimp$vertices)) {
        wimp$vertices[[area_attr]]
      } else {
        rep(NA_character_, nrow(wimp$vertices))
      }

      return(list(
        lpoles = wimp$vertices$left_pole,
        rpoles = wimp$vertices$right_pole,
        wmatrix = wmatrix,
        self = wimp$vertices$self,
        ideal = wimp$vertices$ideal,
        area_vec = area_vec
      ))
    } else if (inherits(wimp, "scn")) {
      return(list(
        lpoles = wimp$constructs[[1]],
        rpoles = wimp$constructs[[2]],
        wmatrix = wimp$weights,
        self = wimp[["self"]][[1]],
        ideal = wimp[["self"]][[2]],
        area_vec = rep(NA_character_, length(wimp$constructs[[1]]))
      ))
    }
    stop("Invalid wimp object type")
  }


  # Create node labels based on vertex values
  .create_vertex_names <- function(vertex_vector, lpoles, rpoles) {
    poles <- paste(lpoles, "-", rpoles)
    sapply(seq_along(vertex_vector), function(n) {
      x <- vertex_vector[n]
      if (x < 0) lpoles[n] else if (x > 0) rpoles[n] else poles[n]
    })
  }

  # Calculate congruency-based colors and groups
  .calculate_congruency <- function(vertex_vector, ideal_vector, color) {
    congruency_vector <- vertex_vector / ideal_vector

    vertex_color <- sapply(congruency_vector, function(x) {
      dplyr::case_when(
        x < 0 && x != -Inf ~ .color_palette(color)[1],  # discrepant
        x > 0 && x !=  Inf ~ .color_palette(color)[2],  # congruent
        x == 0             ~ .color_palette(color)[3],  # undefined
        is.infinite(x)     ~ .color_palette(color)[4],  # dilemmatic
        TRUE               ~ .color_palette(color)[4]   # default to dilemmatic
      )
    })

    vertex_group <- sapply(congruency_vector, function(x) {
      dplyr::case_when(
        x < 0 && x != -Inf ~ "Discrepant",
        x > 0 && x !=  Inf ~ "Congruent",
        x == 0             ~ "Undefined Self",
        is.infinite(x)     ~ "Dilemmatic",
        TRUE               ~ "Dilemmatic"
      )
    })

    list(color = vertex_color, group = vertex_group)
  }

  # Reorient weight matrix based on vertex signs
  .reorient_weight_matrix <- function(wmatrix, vertex_vector) {
    for (n in seq_along(vertex_vector)) {
      if (vertex_vector[n] != 0) {
        direction_value <- vertex_vector[n] / abs(vertex_vector[n])
        wmatrix[, n] <- wmatrix[, n] * direction_value
        wmatrix[n, ] <- wmatrix[n, ] * direction_value
      }
    }
    wmatrix
  }

  # Extract edges from weight matrix
  .extract_edges <- function(wmatrix) {
    n_vertex <- nrow(wmatrix)
    edges_list <- list()

    for (i in seq_len(n_vertex)) {
      for (j in seq_len(n_vertex)) {
        weight <- wmatrix[i, j]
        if (weight != 0) {
          edges_list[[length(edges_list) + 1]] <- data.frame(
            from = i, to = j, weight = weight
          )
        }
      }
    }

    if (length(edges_list) == 0) {
      return(data.frame(from = integer(0), to = integer(0),
                        weight = numeric(0)))
    }

    do.call(rbind, edges_list)
  }

  # Determine edge visual properties based on weight values
  .calculate_edge_properties <- function(weight_vector, color) {
    if (length(weight_vector) == 0) {
      return(list(color = character(0), highlight = character(0), dashes = logical(0)))
    }

    if (color != "grey scale") {
      list(
        color = sapply(weight_vector, function(x) {
          ifelse(x > 0, "rgba(128,128,128,0.5)", "rgba(205,92,92,0.5)")
        }),
        highlight = sapply(weight_vector, function(x) {
          ifelse(x > 0, "rgba(128,128,128,1)", "rgba(205,92,92,1)")
        }),
        dashes = rep(FALSE, length(weight_vector))
      )
    } else {
      list(
        color = rep("rgba(128,128,128,0.5)", length(weight_vector)),
        highlight = rep("rgba(128,128,128,1)", length(weight_vector)),
        dashes = sapply(weight_vector, function(x) {
          ifelse(x > 0, FALSE, TRUE)
        })
      )
    }
  }

  # Detect bidirectional edges for smoothing
  .detect_bidirectional_edges <- function(wmatrix) {
    n_vertex <- nrow(wmatrix)
    edge_curved <- logical()

    for (i in seq_len(n_vertex)) {
      for (j in seq_len(n_vertex)) {
        if (wmatrix[i, j] != 0) {
          # Check if reverse edge exists
          is_bidirectional <- (wmatrix[j, i] != 0)
          edge_curved <- c(edge_curved, is_bidirectional)
        }
      }
    }

    edge_curved
  }

  # ==========================================
  # MAIN FUNCTION LOGIC
  # ==========================================

  # Extract wimp object data
  wimp_data <- .extract_wimp_data(wimp, area_attr)
  lpoles <- wimp_data$lpoles
  rpoles <- wimp_data$rpoles
  wmatrix <- wimp_data$wmatrix
  self <- wimp_data$self
  ideal <- wimp_data$ideal
  area_vec <- wimp_data$area_vec
  poles <- paste(lpoles, "-", rpoles)

  # Validate area attribute
  if (areas && !area_attr %in% names(wimp$vertices) && inherits(wimp, "wimp")) {
    warning("Vertex attribute '", area_attr,
            "' not found. Areas layout will not group by this attribute.")
  }

  # Set default vectors
  if (is.na(vertex_vector[1])) vertex_vector <- self
  vertex_vector <- sapply(vertex_vector, .thr)
  if (is.na(ideal_vector[1])) ideal_vector <- ideal

  # Fix: Adjust weights relative to the displayed pole
  # If we show the left pole (sign < 0), we must flip the incident edges
  # to keep the relationship consistent with the visible label.
  signs <- sign(vertex_vector)
  signs[signs == 0] <- 1 # Default to positive if neutral
  wmatrix <- diag(signs, nrow = length(signs)) %*% wmatrix %*% diag(signs, nrow = length(signs))

  # Create vertex properties
  vertex_name <- .create_vertex_names(vertex_vector, lpoles, rpoles)
  congruency <- .calculate_congruency(vertex_vector, ideal_vector, color)

  # Build vertices data frame
  # "Golden Center" linear formula (balances V1.0 history with modern UI stability)
  # size 20 -> -38 | size 50 -> -71
  vertex_vadjust <- - (33 * abs(vertex_vector) + 38)

  # Color properties for R-side consistency
  border_color <- sapply(congruency$color, function(c) {
    if(c == "#999999") return("#666666")
    if(c == "#FFFF00") return("#CCCC00")
    # Red/Green variants
    if(c == "#F52722") return("#B01C18")
    if(c == "#A5D610") return("#75960B")
    return(c)
  })

  vertex <- data.frame(
    id = as.character(seq_along(vertex_vector)),
    label = vertex_name,
    group = congruency$group,
    category = as.character(area_vec),
    size = 30 * abs(vertex_vector) + 20, 
    shape = "dot",
    title = paste("<p><b>", poles, "</b><br>Self:",
                  round(vertex_vector, 2), "<br>Ideal:",
                  round(ideal_vector, 2), "</p>"),
    color.background = congruency$color,
    color.border = border_color,
    color.highlight.background = congruency$color,
    color.highlight.border = border_color,
    orig_color = congruency$color,
    self_val = as.numeric(vertex_vector),
    ideal_val = as.numeric(ideal_vector),
    raw_size = 30 * abs(vertex_vector) + 20,
    shadow = TRUE,
    hidden = !show,
    font = list(
      face = "Segoe UI", 
      size = 20, 
      color = "#000000", 
      vadjust = vertex_vadjust, 
      strokeWidth = 3, 
      strokeColor = "#ffffff"
    )
  )
  vertex[[area_attr]] <- as.character(area_vec)
  
  # Add other categorical columns from wimp$vertices for tooltips/filtering
  if (inherits(wimp, "wimp")) {
    w_verts <- wimp$vertices
    cat_cols_orig <- names(w_verts)[sapply(w_verts, function(x) is.character(x) || is.factor(x))]
    for (cc in cat_cols_orig) {
      if (!cc %in% names(vertex)) vertex[[cc]] <- as.character(w_verts[[cc]])
    }
  }

  # Internal helper for layout pre-calculation
  .get_all_layouts <- function(wmatrix, vertex, area_attr) {
    ig <- igraph::graph_from_adjacency_matrix(wmatrix, weight = TRUE, mode = "directed")
    
    # igraph layouts
    # Normalize coordinates to a consistent range based on construct count
    .scale <- function(m) {
      if (nrow(m) == 0) return(m)
      # Dynamic multiplier: base 420 + 20 per construct to avoid overcrowding
      mult <- 420 + 20 * nrow(m)
      for (i in 1:2) {
        rng <- range(m[, i])
        span <- rng[2] - rng[1]
        if (span < 1e-9) span <- 1
        m[, i] <- (m[, i] - rng[1]) / span - 0.5
      }
      m * mult
    }
    
    node_ids <- as.character(vertex$id)
    
    layouts <- list(
      "graphopt" = .scale(igraph::layout_with_graphopt(ig)),
      "circle"   = .scale(igraph::layout_in_circle(ig)),
      "mds"      = .scale(igraph::layout_with_mds(ig)),
      "grid"     = .scale(igraph::layout_on_grid(ig)),
      "tree"     = .scale(igraph::layout_as_tree(ig, circular = TRUE))
    )
    
    # Custom Areas layout
    if (areas) {
      v_areas <- .handle_areas_layout(vertex, area_attr)
      layouts[["areas"]] <- .scale(cbind(x = v_areas$x, y = v_areas$y))
    }
    
    # Convert to standard named list for JS with character IDs
    lapply(layouts, function(m) {
      data.frame(id = node_ids, x = as.numeric(m[,1]), y = as.numeric(m[,2]), stringsAsFactors = FALSE)
    })
  }

  # Apply direct relationship hiding if requested
  if (hide_direct && !interactive_options) {
    logical_dilemmatic <- ideal_vector == 0
    wmatrix[wmatrix > 0] <- 0
    wmatrix[logical_dilemmatic, ] <- 0
    wmatrix[, logical_dilemmatic] <- 0
  }

  # Extract all potential edges and their bidirectional status
  all_edges_raw <- .extract_edges(wmatrix)
  all_edge_curved <- .detect_bidirectional_edges(wmatrix)
  
  logical_dilemmatic <- ideal_vector == 0
  
  max_w <- (if (nrow(all_edges_raw) > 0) max(abs(all_edges_raw$weight)) else 1) + 0.01

  edges_raw <- all_edges_raw
  edge_curved <- all_edge_curved

  if (nrow(edges_raw) == 0) {
    edges <- data.frame(
      from = integer(0), to = integer(0), width = numeric(0),
      arrows = character(0), dashes = logical(0),
      color.color = character(0), color.highlight = character(0), title = numeric(0), weight = numeric(0),
      hidden = logical(0), is_dilemmatic = logical(0),
      stringsAsFactors = FALSE
    )
  } else {
    edge_props <- .calculate_edge_properties(edges_raw$weight, color)
    
    is_dilemmatic_edge <- logical_dilemmatic[edges_raw$from] | logical_dilemmatic[edges_raw$to]

    edges <- data.frame(
      from = as.character(edges_raw$from),
      to = as.character(edges_raw$to),
      width = 2 * abs(edges_raw$weight),
      arrows = "to",
      dashes = edge_props$dashes,
      color.color = edge_props$color,
      color.highlight = edge_props$highlight,
      title = round(edges_raw$weight, 2),
      weight = edges_raw$weight,
      orig_dashes = edge_props$dashes,
      orig_color = edge_props$color,
      orig_highlight = edge_props$highlight,
      is_direct = edges_raw$weight > 0,
      is_dilemmatic = is_dilemmatic_edge,
      hidden = if (interactive_options) (abs(edges_raw$weight) < min_weight | (hide_direct & (edges_raw$weight > 0 | is_dilemmatic_edge))) else FALSE,
      stringsAsFactors = FALSE
    )
  }

  # ==========================================
  # LAYOUT PROCESSING
  # ==========================================

  # Convert layout names to igraph format
  .convert_layout_name <- function(layout) {
    layout_mapping <- list(
      "graphopt" = "layout_with_graphopt",
      "circle"   = "layout_in_circle",
      "tree"     = "layout_as_tree",
      "mds"      = "layout_with_mds",
      "grid"     = "layout_on_grid",
      "rtcircle" = "layout_in_circle" # rtcircle uses layout_in_circle with circular=TRUE
    )
    if (layout %in% names(layout_mapping)) {
      layout_mapping[[layout]]
    } else {
      layout
    }
  }

  # Handle areas-specific layout
  .handle_areas_layout <- function(vertex, area_attr) {
    cat2 <- vertex[[area_attr]]
    cat2[is.na(cat2) | cat2 == ""] <- "Uncategorized"
    vertex[[area_attr]] <- cat2
    cats <- unique(cat2)
    k <- length(cats)

    # Grid layout constants
    ncol <- ceiling(sqrt(k))
    nrow <- ceiling(k / ncol)
    cell <- 400

    # Calculate grid centers
    centers <- lapply(seq_len(k), function(i) {
      r <- floor((i - 1) / ncol)
      c <- (i - 1) %% ncol
      c(cx = (c - (ncol - 1) / 2) * cell,
        cy = ((nrow - 1) / 2 - r) * cell)
    })
    centers <- do.call(rbind, centers)
    rownames(centers) <- cats

    # Position nodes within each category
    pos_list <- lapply(cats, function(cat) {
      ids <- vertex$id[vertex[[area_attr]] == cat & !vertex$hidden]
      if (length(ids) == 0L) return(NULL)
      cx <- centers[cat, "cx"]
      cy <- centers[cat, "cy"]
      .place_in_cell(ids, cx, cy, cell)
    })

    pos_df <- do.call(rbind, pos_list)
    vertex$x <- NA_real_
    vertex$y <- NA_real_

    if (!is.null(pos_df)) {
      vertex$x[match(pos_df$id, vertex$id)] <- pos_df$x
      vertex$y[match(pos_df$id, vertex$id)] <- pos_df$y
    }

    # Handle unplaced nodes
    missing <- is.na(vertex$x)
    if (any(missing)) {
      set.seed(33)
      vertex$x[missing] <- rnorm(sum(missing), 0, cell * 0.05)
      vertex$y[missing] <- rnorm(sum(missing), 0, cell * 0.05)
    }

    vertex
  }

  # ==========================================
  # NETWORK RENDERING
  # ==========================================

  # Determine layout type and process accordingly
  use_areas_layout <- identical(layout, "areas")
  if (use_areas_layout) {
    vertex <- .handle_areas_layout(vertex, area_attr)
  }

  layout_name <- .convert_layout_name(layout)

  # Create base network
  if (use_areas_layout) {
    # Areas layout with fixed positions
    g <- visNetwork(vertex, edges, height = height, width = width) %>%
      visOptions(manipulation = list(enabled = FALSE),
                 highlightNearest = list(enabled = TRUE, degree = 0,
                                         labelOnly = TRUE)) %>%
      visInteraction(navigationButtons = FALSE, multiselect = TRUE, selectConnectedEdges = FALSE) %>%
      visPhysics(enabled = FALSE) %>%
      visEdges(smooth = list(enabled = TRUE, type = "curvedCW", roundness = 0.15)) %>% 
      visNodes(font = list(align = "center", multi = TRUE, vadjust = 0))
  } else if (layout == "rtcircle") {
    # Circular tree layout
    g <- visNetwork(vertex, edges, height = height, width = width) %>%
      visIgraphLayout(layout = "layout_as_tree", circular = TRUE) %>%
      visOptions(manipulation = list(enabled = FALSE),
                 highlightNearest = list(enabled = TRUE, degree = 0,
                                         labelOnly = TRUE)) %>%
      visInteraction(navigationButtons = FALSE, multiselect = TRUE, selectConnectedEdges = FALSE) %>%
      visEdges(smooth = list(enabled = TRUE, type = "curvedCW", roundness = 0.15)) %>% 
      visNodes(font = list(align = "center", multi = TRUE, vadjust = 0))
  } else {
    # Standard igraph layouts
    g <- visNetwork(vertex, edges, height = height, width = width) %>%
      visIgraphLayout(layout = layout_name, randomSeed = 33) %>%
      visOptions(manipulation = list(enabled = FALSE),
                 highlightNearest = list(enabled = TRUE, degree = 0,
                                         labelOnly = TRUE)) %>%
      visInteraction(navigationButtons = FALSE, multiselect = TRUE, selectConnectedEdges = FALSE) %>%
      visEdges(smooth = list(enabled = TRUE, type = "curvedCW", roundness = 0.15)) %>% 
      visNodes(font = list(align = "center", multi = TRUE, vadjust = 0))
  }

  # ==========================================
  # VIEWPORT AND INTERACTION SETUP
  # ==========================================

  # Configure auto-centering behavior
  if (!use_areas_layout || !areas) {
    g <- g %>% visEvents(
      stabilized = htmlwidgets::JS(
        "function(){ if(!this._autoFitDone){ try{ this.fit(); }catch(e){} ",
        "this._autoFitDone=true; } }"
      ),
      afterDrawing = htmlwidgets::JS(
        "function(ctx){ if(!this._autoFitDone){ try{ this.fit(); }catch(e){} ",
        "this._autoFitDone=true; } }"
      )
    )
  } else {
    # For areas layout, fit once after stabilization
    g <- g %>% visEvents(
      stabilized = htmlwidgets::JS(
        "function(){ if(!this._autoFitDone){ try{ this.fit({animation: ",
        "false}); }catch(e){} this._autoFitDone=true; } }"
      )
    )
  }

  # Ensure category attribute is preserved in widget nodes
  if ("category" %in% names(vertex) && !is.null(g$x$nodes)) {
    if (NROW(g$x$nodes) == nrow(vertex)) {
      g$x$nodes$category <- vertex$category
    }
  }

  # ==========================================
  # AREA OVERLAYS (CONVEX HULLS)
  # ==========================================

  if (areas) {
    uniq_cats <- unique(stats::na.omit(vertex$category[vertex$category != ""]))
    if (length(uniq_cats) > 0) {
      # Define area colors
      if (missing(area_color) || is.null(area_color) ||
            all(is.na(area_color))) {
        base_cols <- c("#6EA8FE", "#72D6A0", "#F7B267", "#D985B9",
                       "#8FD3FE", "#B5E48C", "#FFD166", "#CDB4DB",
                       "#A0C4FF", "#FFAFCC")
      } else {
        base_cols <- area_color
      }

      # Calculate color vectors and transparency
      fill_vec   <- rep(base_cols, length.out = length(uniq_cats))
      stroke_vec <- fill_vec

      # Convert hex colors to rgba with transparency
      rgb_to_rgba <- function(hex, a) {
        hex <- gsub("#", "", hex)
        r <- strtoi(substr(hex, 1, 2), 16)
        g <- strtoi(substr(hex, 3, 4), 16)
        b <- strtoi(substr(hex, 5, 6), 16)
        sprintf("rgba(%d,%d,%d,%.2f)", r, g, b, a)
      }

      fills   <- vapply(fill_vec,   rgb_to_rgba, character(1), a = 0.15)
      strokes <- vapply(stroke_vec, rgb_to_rgba, character(1), a = 0.65)

      # Generate JavaScript for convex hulls
      js_hulls <- .create_js_hulls(pad_side, rounding, uniq_cats, fills,
                                   strokes)
      g <- g %>% visEvents(afterDrawing = htmlwidgets::JS(js_hulls))
    }
  }

  # ==========================================
  # DYNAMIC CONTROLS (Unified Graph Options Panel)
  # ==========================================
  
  if (interactive_options) {
    js_panel <- "
    function(el, x) {
      el.style.height = '100%';
      el.style.minHeight = '0px';
      var network = this.network;
      var container = el;
      el.style.position = 'relative';

      // Inject Chart.js if not present
      if (!window.Chart) {
        var script = document.createElement('script');
        script.src = 'https://cdn.jsdelivr.net/npm/chart.js';
        document.head.appendChild(script);
      }
      
      var currentDistMult = 1.0;
      var currentSizeMult = 1.0;
      var currentTextSize = 20;
      var viewMode = 'graph';
      var pcsdChart = null;

      // --- Utilities ---
      var getPaletteColor = function(v, i, scheme) {
        var p = x.color_palette_js[scheme] || x.color_palette_js['red/green'];
        if (i === 0) return p[3]; 
        if (v === 0) return p[2];
        return (Math.sign(v) === Math.sign(i)) ? p[1] : p[0];
      };

      var darkenColor = function(hex, percent) {
        hex = hex.replace('#', '');
        var r = parseInt(hex.substring(0, 2), 16),
            g = parseInt(hex.substring(2, 4), 16),
            b = parseInt(hex.substring(4, 6), 16);
        r = Math.floor(r * (1 - percent));
        g = Math.floor(g * (1 - percent));
        b = Math.floor(b * (1 - percent));
        return '#' + ((1 << 24) + (r << 16) + (g << 8) + b).toString(16).slice(1);
      };

      var refreshNodes = function(simVals) {
        var scheme = visContent.querySelector('#palette_sel').value;
        var currentIdx = (network._simCurrentI !== undefined) ? network._simCurrentI : 0;
        var vals = simVals || (network._simHistory && network._simHistory[currentIdx]);
        var nodesDS = network.body.data.nodes;
        
        var updates = nodesDS.get().map(function(node) {
          var idx = (parseInt(node.id) - 1);
          if (isNaN(idx)) return node;
          var val = (vals && vals[idx] !== undefined) ? vals[idx] : (node.self_val || 0);
          var ideal = (x.sim_data && x.sim_data.initial_ideal) ? x.sim_data.initial_ideal[idx] : node.ideal_val;
          var c = getPaletteColor(val, ideal, scheme);
          var label = node.label;
          if (x.sim_data && network._simHistory && x.sim_data.lpoles && x.sim_data.rpoles) {
            var lp = x.sim_data.lpoles[idx] || 'L';
            var rp = x.sim_data.rpoles[idx] || 'R';
            label = (val < -0.1) ? lp : (val > 0.1 ? rp : lp + ' - ' + rp);
          }
          var baseSize = (x.sim_data && network._simHistory) ? (20 + (30 * Math.abs(val))) : (node.raw_size || 20);
          var finalSize = baseSize * currentSizeMult;
          var vadjust = - (finalSize * 1.1 + currentTextSize * 0.8);
          return {
            id: node.id,
            color: { background: c, border: darkenColor(c, 0.4), highlight: { background: c, border: darkenColor(c, 0.4) } },
            size: finalSize, label: label, shape: 'dot',
            font: { vadjust: vadjust, size: currentTextSize, face: 'Segoe UI', color: '#000000', strokeWidth: 3, strokeColor: '#ffffff' }
          };
        });
        nodesDS.update(updates);
        
        var edgesDS = network.body.data.edges;
        var scheme2 = visContent.querySelector('#palette_sel').value;
        var opSlider = visContent.querySelector('#opacity_slider');
        var op = opSlider ? parseFloat(opSlider.value) : 0.5;
        
        var edgeUpdates = edgesDS.get().map(function(edge) {
          var c1 = edge.orig_color ? edge.orig_color.replace(/[\\d.]+\\)$/, op + ')') : '';
          var c2 = (scheme2 === 'grey scale') ? 'rgba(153,153,153,'+op+')' : c1;
          var hl = (scheme2 === 'grey scale') ? 'rgba(153,153,153,1)' : edge.orig_highlight;
          var dsh = (scheme2 === 'grey scale') ? (edge.weight < 0) : edge.orig_dashes;
          return {id: edge.id, color: {color: c2, highlight: hl}, dashes: dsh};
        });
        edgesDS.update(edgeUpdates);
      };

      // ── Semantic Edge flow flash: arrows showing activation quality ───────
      var flashEdgeFlow = function(fromIdx, toIdx) {
        network._activePulses = [];
        if(fromIdx === 0) return; // Step 0->1 is static setup
        
        // Ripple Effect: Pulse based on the change that happened in the PREVIOUS step (n-1 to n)
        if(!network._simHistory || !network._simHistory[fromIdx-1] || !network._simHistory[fromIdx] || !x.sim_data) return;
        var prevHistory = network._simHistory[fromIdx - 1]; 
        var currHistory = network._simHistory[fromIdx];     
        var deltas = currHistory.map(function(v, i) { return v - prevHistory[i]; });
        
        var weights  = x.sim_data.weights;
        var ideals   = x.sim_data.initial_ideal || [];
        var THRESHOLD = 0.04;
        var edgesDS  = network.body.data.edges;
        
        var flowUpdates = edgesDS.get().map(function(edge) {
          var srcIdx = parseInt(edge.from) - 1;
          var dstIdx = parseInt(edge.to) - 1;
          if(isNaN(srcIdx) || isNaN(dstIdx)) return {id: edge.id};
          
          var flow = 0;
          if(weights && weights[srcIdx]) {
             if (Array.isArray(weights[srcIdx])) flow = weights[srcIdx][dstIdx] * deltas[srcIdx];
             else flow = weights[srcIdx * ideals.length + dstIdx] * deltas[srcIdx]; 
          }
          
          if (Math.abs(flow) > THRESHOLD) {
            var idealVal = ideals[dstIdx] || 0;
            var type, color;
            if (idealVal !== 0) {
              var isToward = (Math.sign(flow) === Math.sign(idealVal));
              type  = isToward ? 'up' : 'down';
              color = isToward ? '#4CAF50' : '#E53935';
            } else {
              type  = (flow > 0) ? 'right' : 'left';
              color = '#FBC02D';
            }
            network._activePulses.push({ from: edge.from, to: edge.to, color: color, type: type, t: 0 });
            return {id: edge.id, color: {color: color, highlight: color, hover: color}};
          }
          return {id: edge.id, color: edge.orig_color || '#cccccc'};
        });
        edgesDS.update(flowUpdates);
      };

      // ── Arrow Pulse Renderer: Drawing glowing arrows ────────────────────
      network.on(\"afterDrawing\", function(ctx) {
        if (!network._activePulses || network._activePulses.length === 0) return;
        var positions = network.getPositions();
        network._activePulses.forEach(function(p) {
          var start = positions[p.from];
          var end   = positions[p.to];
          if (!start || !end) return;
          var posX = start.x + (end.x - start.x) * p.t;
          var posY = start.y + (end.y - start.y) * p.t;
          
          ctx.save();
          ctx.translate(posX, posY);
          ctx.shadowBlur = 20; // Thicker glow
          ctx.shadowColor = p.color;
          ctx.fillStyle = p.color;
          
          var sz = 9; // Slightly larger
          var drawArrowShape = function(s) {
            ctx.beginPath();
            var sw = s/1.8; // slightly thicker stem
            if(p.type === 'up') {
              ctx.moveTo(-sw/2, s); ctx.lineTo(sw/2, s); ctx.lineTo(sw/2, 0); ctx.lineTo(s, 0); ctx.lineTo(0, -s); ctx.lineTo(-s, 0); ctx.lineTo(-sw/2, 0);
            } else if(p.type === 'down') {
              ctx.moveTo(-sw/2, -s); ctx.lineTo(sw/2, -s); ctx.lineTo(sw/2, 0); ctx.lineTo(s, 0); ctx.lineTo(0, s); ctx.lineTo(-s, 0); ctx.lineTo(-sw/2, 0);
            } else if(p.type === 'right') {
              ctx.moveTo(-s, -sw/2); ctx.lineTo(-s, sw/2); ctx.lineTo(0, sw/2); ctx.lineTo(0, s); ctx.lineTo(s, 0); ctx.lineTo(0, -s); ctx.lineTo(0, -sw/2);
            } else { // left
              ctx.moveTo(s, -sw/2); ctx.lineTo(s, sw/2); ctx.lineTo(0, sw/2); ctx.lineTo(0, s); ctx.lineTo(-s, 0); ctx.lineTo(0, -s); ctx.lineTo(0, -sw/2);
            }
            ctx.closePath();
            ctx.fill();
          };
          
          drawArrowShape(sz);
          ctx.shadowBlur = 0;
          ctx.fillStyle = '#ffffff';
          drawArrowShape(sz/2); // inner core relative to size
          ctx.restore();
        });
      });

      // ── Dynamic Area Drawing (Convex Hulls) ────────────────────────────────
      var drawAreasOnCanvas = function(ctx) {
        if(currentAreaAttr === 'None') return;
        
        var PAD = 50; 
        var ROUND = 10;
        var LINE_W = 2;
        var MAX_OFFSET_DIST = 150;
        
        var FILL = ['rgba(110,168,254,0.15)', 'rgba(114,214,160,0.15)', 'rgba(247,178,103,0.15)', 'rgba(217,133,185,0.15)', 'rgba(143,211,254,0.15)', 'rgba(181,228,140,0.15)', 'rgba(255,209,102,0.15)', 'rgba(205,180,219,0.15)', 'rgba(160,196,255,0.15)', 'rgba(255,175,204,0.15)'];
        var STROK = ['rgba(110,168,254,0.65)', 'rgba(114,214,160,0.65)', 'rgba(247,178,103,0.65)', 'rgba(217,133,185,0.65)', 'rgba(143,211,254,0.65)', 'rgba(181,228,140,0.65)', 'rgba(255,209,102,0.65)', 'rgba(205,180,219,0.65)', 'rgba(160,196,255,0.65)', 'rgba(255,175,204,0.65)'];
        
        function isCCW(pts){var s=0;for(var i=0;i<pts.length;i++){var a=pts[i],b=pts[(i+1)%pts.length];s+=a.x*b.y-a.y*b.x;}return s>0;}
        function distance(p1,p2){return Math.hypot(p2.x-p1.x,p2.y-p1.y);}
        function lineIntersect(p1,d1,p2,d2){var det=d1.x*d2.y-d1.y*d2.x;if(Math.abs(det)<1e-9) return null;var t=((p2.x-p1.x)*d2.y-(p2.y-p1.y)*d2.x)/det;return {x:p1.x+t*d1.x,y:p1.y+t*d1.y};}
        function convexHull(points){
          if(points.length<=1) return points.slice();
          var pts=points.slice().sort(function(a,b){return a.x!==b.x?a.x-b.x:a.y-b.y;});
          function cross(o,a,b){return (a.x-o.x)*(b.y-o.y)-(a.y-o.y)*(b.x-o.x);}
          var lower=[],upper=[];
          for(var i=0;i<pts.length;i++){while(lower.length>=2 && cross(lower[lower.length-2], lower[lower.length-1], pts[i])<=0) lower.pop(); lower.push(pts[i]);}
          for(var j=pts.length-1;j>=0;j--){while(upper.length>=2 && cross(upper[upper.length-2], upper[upper.length-1], pts[j])<=0) upper.pop(); upper.push(pts[j]);}
          upper.pop(); lower.pop(); return lower.concat(upper);
        }
        function offsetConvex(pts, d){
          if(pts.length<3) return null;
          var H=pts.slice(); if(!isCCW(H)) H.reverse();
          var N=H.length, out=new Array(N);
          for(var i=0;i<N;i++){
            var p0=H[(i-1+N)%N], p1=H[i], p2=H[(i+1)%N];
            var e0={x:p1.x-p0.x,y:p1.y-p0.y}, e1={x:p2.x-p1.x,y:p2.y-p1.y};
            var l0=Math.hypot(e0.x,e0.y)||1, l1=Math.hypot(e1.x,e1.y)||1; e0.x/=l0; e0.y/=l0; e1.x/=l1; e1.y/=l1;
            var n0={x:e0.y,y:-e0.x}, n1={x:e1.y,y:-e1.x};
            var dot=e0.x*e1.x+e0.y*e1.y; var angle=Math.acos(Math.max(-1,Math.min(1,-dot)));
            var isSharpCorner=angle<Math.PI/3;
            var pB={x:p1.x+n0.x*d,y:p1.y+n0.y*d}, pC={x:p1.x+n1.x*d,y:p1.y+n1.y*d};
            var q=lineIntersect(pB,e0,pC,e1);
            if(q && distance(p1,q)<=MAX_OFFSET_DIST && !isSharpCorner){
              out[i]=q;
            } else {
              var avgNormal={x:(n0.x+n1.x)*0.5, y:(n0.y+n1.y)*0.5};
              var len=Math.hypot(avgNormal.x,avgNormal.y)||1;
              avgNormal.x/=len; avgNormal.y/=len;
              var safeDist=isSharpCorner ? d*1.5 : d;
              out[i]={x:p1.x+avgNormal.x*safeDist, y:p1.y+avgNormal.y*safeDist};
            }
          }
          return out;
        }
        function draw(ctx, pts, fill, stroke){
          if(!pts || pts.length<3) return;
          ctx.save(); 
          ctx.fillStyle=fill; ctx.strokeStyle=stroke; ctx.lineWidth=LINE_W;
          if(ROUND<=0){
            ctx.beginPath(); ctx.moveTo(pts[0].x,pts[0].y);
            for(var i=1;i<pts.length;i++) ctx.lineTo(pts[i].x,pts[i].y);
            ctx.closePath(); ctx.fill(); ctx.stroke();
          }else{
            ctx.lineJoin='round'; ctx.lineCap='round'; ctx.miterLimit=4;
            ctx.beginPath(); var n=pts.length;
            for(var i=0;i<n;i++){
              var p0=pts[(i-1+n)%n], p1=pts[i], p2=pts[(i+1)%n];
              var v1x=p1.x-p0.x,v1y=p1.y-p0.y,v2x=p2.x-p1.x,v2y=p2.y-p1.y;
              var l1=Math.hypot(v1x,v1y)||1,l2=Math.hypot(v2x,v2y)||1;
              var dot=(v1x/l1)*(-v2x/l2)+(v1y/l1)*(-v2y/l2);
              var angle=Math.acos(Math.max(-1,Math.min(1,dot)));
              var angleRatio=Math.max(0.3,Math.sin(angle*0.5));
              var adaptiveRound=ROUND*angleRatio;
              var rr=Math.min(adaptiveRound,0.4*l1,0.4*l2);
              v1x/=l1; v1y/=l1; v2x/=l2; v2y/=l2;
              var p1_in={x:p1.x-v1x*rr,y:p1.y-v1y*rr}, p1_out={x:p1.x+v2x*rr,y:p1.y+v2y*rr};
              if(i===0) ctx.moveTo(p1_in.x,p1_in.y); else ctx.lineTo(p1_in.x,p1_in.y);
              if(angle<Math.PI/4){
                var ctrl={x:(p1_in.x+p1.x+p1_out.x)/3, y:(p1_in.y+p1.y+p1_out.y)/3};
                ctx.quadraticCurveTo(ctrl.x,ctrl.y,p1_out.x,p1_out.y);
              } else {
                ctx.arcTo(p1.x,p1.y,p1_out.x,p1_out.y,rr);
              }
            }
            ctx.closePath(); ctx.fill(); ctx.stroke();
          }
          ctx.restore();
        }
        function drawLabelAbove(ctx, pts, text, stroke){
          var xs = pts.map(function(p){ return p.x; }), ys = pts.map(function(p){ return p.y; });
          var minX = Math.min.apply(null,xs), maxX = Math.max.apply(null,xs);
          var minY = Math.min.apply(null,ys);
          var cx = (minX + maxX) / 2, y = minY - 12;
          ctx.save(); ctx.font='bold 14px sans-serif'; ctx.textAlign='center'; ctx.textBaseline='bottom';
          ctx.strokeStyle = stroke; ctx.lineWidth = 4; ctx.strokeText(text, cx, y);
          ctx.fillStyle = 'rgba(0,0,0,0.80)'; ctx.fillText(text, cx, y); ctx.restore();
        }

        var nodes=network.body.data.nodes.get(); 
        var bycat={};
        nodes.forEach(function(n){
          if(n.hidden) return;
          if(n[currentAreaAttr]==null) return;
          var c=String(n[currentAreaAttr]);
          if(c===''||c==='NA') return;
          (bycat[c]||(bycat[c]=[])).push(n.id);
        });

        var CAT = Object.keys(bycat);
        for(var i=0;i<CAT.length;i++){
          var cat = CAT[i], ids = bycat[cat];
          if(!ids || ids.length===0) continue;
          var pos = network.getPositions(ids);
          var pts = ids.map(function(id){ return {x:pos[id].x, y:pos[id].y}; });

          var outline;
          if(pts.length<=2){
            var xs=pts.map(function(p){ return p.x; }), ys=pts.map(function(p){ return p.y; });
            var minX=Math.min.apply(null,xs), maxX=Math.max.apply(null,xs);
            var minY=Math.min.apply(null,ys), maxY=Math.max.apply(null,ys);
            var width=maxX-minX, height=maxY-minY;
            var adaptivePad=Math.max(PAD, Math.max(width,height)*0.3+20);
            outline=[
              {x:minX-adaptivePad,y:minY-adaptivePad},
              {x:maxX+adaptivePad,y:minY-adaptivePad},
              {x:maxX+adaptivePad,y:maxY+adaptivePad},
              {x:minX-adaptivePad,y:maxY+adaptivePad}
            ];
          } else {
            var xs=pts.map(function(p){ return p.x; }), ys=pts.map(function(p){ return p.y; });
            var avgX=xs.reduce(function(a,b){return a+b;})/xs.length, avgY=ys.reduce(function(a,b){return a+b;})/ys.length;
            var avgDist=pts.reduce(function(sum,p){return sum+distance({x:avgX,y:avgY},p);},0)/pts.length;
            var adaptivePad=Math.max(PAD, avgDist*0.15+15);
            var hull=convexHull(pts); 
            outline=offsetConvex(hull, adaptivePad);
            if(!outline){
              var minX=Math.min.apply(null,xs)-adaptivePad, maxX=Math.max.apply(null,xs)+adaptivePad;
              var minY=Math.min.apply(null,ys)-adaptivePad, maxY=Math.max.apply(null,ys)+adaptivePad;
              outline=[{x:minX,y:minY},{x:maxX,y:minY},{x:maxX,y:maxY},{x:minX,y:maxY}];
            }
          }

          var colorIdx = i % FILL.length;
          draw(ctx, outline, FILL[colorIdx], STROK[colorIdx]);
          drawLabelAbove(ctx, outline, String(cat), STROK[colorIdx]);
        }
      };

      network.on('beforeDrawing', function(ctx) {
        drawAreasOnCanvas(ctx);
      });

      var _tweenRAF = null;
      var tweenToIteration = function(fromIdx, toIdx, durationMs) {
        if(_tweenRAF) { cancelAnimationFrame(_tweenRAF); _tweenRAF = null; }
        if(!network._simHistory || !network._simHistory[fromIdx] || !network._simHistory[toIdx]) {
          network._simCurrentI = toIdx;
          refreshNodes();
          return;
        }
        var prevVals = network._simHistory[fromIdx].slice();
        var nextVals = network._simHistory[toIdx];
        var nodesDS  = network.body.data.nodes;
        var allNodes = nodesDS.get();
        var scheme   = visContent.querySelector('#palette_sel').value;
        var t0 = null;
        flashEdgeFlow(fromIdx, toIdx);
        var step = function(ts) {
          if(!t0) t0 = ts;
          var t = Math.min((ts - t0) / durationMs, 1);
          var isInit = fromIdx === 0;
          
          // Phase 1: Arrows travel (only if not init Step 0->1)
          var pulseDur = isInit ? 0 : 0.7;
          var pulseT = (pulseDur === 0) ? 0 : Math.min(t / pulseDur, 1);
          if(network._activePulses) {
            network._activePulses.forEach(function(p) { p.t = pulseT; });
          }

          // Phase 2: Nodes update (starts after arrows if not init)
          var nodeStart = isInit ? 0 : 0.7;
          var nodeT = t < nodeStart ? 0 : (t - nodeStart) / (1 - nodeStart);
          var easeUpdate = nodeT < 0.5 ? 4*nodeT*nodeT*nodeT : 1 - Math.pow(-2*nodeT+2, 3)/2;
          
          var interpVals = nextVals.map(function(nv, i) { return prevVals[i] + easeUpdate * (nv - prevVals[i]); });
          var updates = allNodes.map(function(node) {
            var idx   = parseInt(node.id) - 1;
            if (isNaN(idx)) return node;
            var val   = (interpVals && interpVals[idx] !== undefined) ? interpVals[idx] : (node.self_val || 0);
            var ideal = (x.sim_data && x.sim_data.initial_ideal) ? x.sim_data.initial_ideal[idx] : node.ideal_val;
            var c     = getPaletteColor(val, ideal, scheme);
            var baseSize = 20 + 30 * Math.abs(val);
            var finalSize = baseSize * currentSizeMult;
            var vadjust   = -(finalSize * 1.1 + currentTextSize * 0.8);
            var label = node.label;
            if (x.sim_data && x.sim_data.lpoles && x.sim_data.rpoles) {
               var lp = x.sim_data.lpoles[idx] || 'L';
               var rp = x.sim_data.rpoles[idx] || 'R';
               label = (val < -0.1) ? lp : (val > 0.1 ? rp : lp + ' - ' + rp);
            }
            return {
              id: node.id,
              color: { background: c, border: darkenColor(c, 0.4), highlight: { background: c, border: darkenColor(c, 0.6) } },
              size: finalSize, label: label, shape: 'dot',
              font: { vadjust: vadjust, size: currentTextSize, face: 'Segoe UI', color: '#000000', strokeWidth: 3, strokeColor: '#ffffff' }
            };
          });
          nodesDS.update(updates); // This triggers a redraw which calls afterDrawing
          if(t < 1) {
            _tweenRAF = requestAnimationFrame(step);
          } else {
            _tweenRAF = null;
            network._simCurrentI = toIdx;
            network._activePulses = []; // Clear pulses
            refreshNodes();
          }
        };
        _tweenRAF = requestAnimationFrame(step);
      };


      // --- Export Logic ---
      var exportPNG = function() {
        var canvas = container.getElementsByTagName('canvas')[0];
        if(!canvas) return;
        var link = document.createElement('a');
        link.download = 'WIMP_EXPORT_NAME_Digraph.png';
        link.href = canvas.toDataURL('image/png', 1.0);
        link.click();
      };

      var isSidebarMode = !!document.getElementById('wsim_sidebar');
      // --- Helper to create panels ---
      var createPanel = function(id, title, positionStyles) {
        var p = document.createElement('div');
        p.id = id;
        var useSidebar = (isSidebarMode && id === 'sim_settings_panel');
        var sidebar = document.getElementById('wsim_sidebar');
        if (useSidebar) {
          Object.assign(p.style, {
            backgroundColor: '#ffffff',
            padding: '15px', borderRadius: '8px', boxShadow: '0 2px 10px rgba(0,0,0,0.05)',
            border: '1px solid #ddd', fontFamily: 'Segoe UI, Tahoma, sans-serif', fontSize: '12px',
            marginBottom: '0', transition: 'all 0.3s ease', width: '100%',
            display: 'flex', flexDirection: 'column', flex: '1', boxSizing: 'border-box', minHeight: '0'
          });
          p.style.maxHeight = 'none';
        } else {
          Object.assign(p.style, {
            position: 'absolute', zIndex: '1000', backgroundColor: 'rgba(255, 255, 255, 0.95)',
            padding: '10px', borderRadius: '8px', boxShadow: '0 2px 15px rgba(0,0,0,0.15)',
            border: '1px solid #ddd', fontFamily: 'Segoe UI, Tahoma, sans-serif', fontSize: '12px',
            width: '90%', maxWidth: '220px', boxSizing: 'border-box', maxHeight: '40px', overflowY: 'hidden', transition: 'all 0.3s ease'
          }, positionStyles);
        }

        var header = document.createElement('div');
        header.style.display = 'flex';
        header.style.justifyContent = 'space-between';
        header.style.alignItems = 'center';
        header.style.cursor = 'pointer';
        
        var titleTxt = document.createElement('b');
        titleTxt.style.color = '#2c3e50';
        titleTxt.innerText = title;
        header.appendChild(titleTxt);

        var toggleIcon = document.createElement('span');
        if (!useSidebar) {
            toggleIcon.className = 'toggle-icon';
            toggleIcon.innerText = '+';
            header.appendChild(toggleIcon);
        } else {
            header.style.cursor = 'default';
        }
        p.appendChild(header);

        var content = document.createElement('div');
        content.style.display = 'none';
        content.style.marginTop = '10px';
        content.style.borderTop = '1px solid #eee';
        content.style.paddingTop = '10px';
        p.appendChild(content);

        if (useSidebar) {
          content.style.display = 'flex';
          content.style.flexDirection = 'column';
          content.style.flex = '1';
          content.style.overflow = 'hidden';
          content.style.minHeight = '0';
        }

        header.onclick = function() {
          if (useSidebar) return; // Do not collapse in sidebar
          var isHidden = content.style.display === 'none';
          content.style.display = isHidden ? 'block' : 'none';
          toggleIcon.innerText = isHidden ? '−' : '+';
          p.style.maxHeight = isHidden ? '90%' : '40px';
          p.style.overflowY = isHidden ? 'auto' : 'hidden';
        };

        if (useSidebar) {
          sidebar.appendChild(p);
        } else {
          container.appendChild(p);
        }
        return content;
      };

      // --- Settings Modal Initialization ---
      var settingsModal = document.createElement('div');
      settingsModal.id = 'settings_modal';
      Object.assign(settingsModal.style, {
        position: 'absolute', top: '10px', right: '10px',
        width: '90%', maxWidth: '300px', maxHeight: '80vh', overflowY: 'auto', boxSizing: 'border-box',
        backgroundColor: '#fff', zIndex: '2000', padding: '20px', borderRadius: '8px', 
        boxShadow: '0 4px 20px rgba(0,0,0,0.2)', border: '1px solid #eaeaea', display: 'none', 
        fontFamily: 'Inter, Roboto, sans-serif'
      });
      
      var settingsModalHeader = '<div style=\"display:flex; justify-content:space-between; align-items:center; border-bottom:1px solid #eaeaea; padding-bottom:10px; margin-bottom:15px;\">' +
                                '<h3 style=\"margin:0; color:#444; font-size:14px;\">' + x.dict.vis_options + '</h3>' +
                                '<span id=\"close_settings_modal\" style=\"cursor:pointer; font-size:20px; font-weight:bold; color:#888; line-height:1;\">&times;</span>' +
                                '</div>';
                                
      var settingsModalBody = document.createElement('div');
      settingsModal.innerHTML = settingsModalHeader;
      settingsModal.appendChild(settingsModalBody);
      container.appendChild(settingsModal);
      
      settingsModal.querySelector('#close_settings_modal').onclick = function() { settingsModal.style.display = 'none'; };

      var visContent = document.createElement('div');
      
      // --- Area Selector Panel ---
      var currentAreaAttr = 'None';
      if (x.cat_cols && x.cat_cols.length > 0) {
        var areaPanelContent = createPanel('area_panel', x.dict.areas, {top: '10px', left: '10px'});
        
        var areaHTML = '<div style=\"margin-bottom:10px;\">' +
                       '<label style=\"display:block; margin-bottom:5px; font-weight:bold; color:#444;\">' + x.dict.group_by + '</label>' +
                       '<select id=\"area_sel\" style=\"width:100%; padding:4px; border-radius:4px; margin-bottom:8px;\">' +
                       '<option value=\"None\">' + x.dict.none + '</option>' +
                       x.cat_cols.map(function(c) { return '<option value=\"' + c + '\">' + c.charAt(0).toUpperCase() + c.slice(1) + '</option>'; }).join('') +
                       '</select>' +
                       '<button id=\"btn_area_layout\" style=\"width:100%; padding:6px; background:#f4f9ef; border:1px solid #8cc63f; color:#5c8822; border-radius:4px; cursor:pointer; font-size:11px; font-weight:bold; display:none; align-items:center; justify-content:center; gap:4px;\">' + x.dict.cluster_areas + '</button>' +
                       '</div>';
        areaPanelContent.innerHTML = areaHTML;
        
        areaPanelContent.querySelector('#area_sel').onchange = function() {
          currentAreaAttr = this.value;
          areaPanelContent.querySelector('#btn_area_layout').style.display = (currentAreaAttr === 'None') ? 'none' : 'flex';
          network.redraw();
        };
        
        areaPanelContent.querySelector('#btn_area_layout').onclick = function() {
          if (currentAreaAttr === 'None') return;
          var nodes = network.body.data.nodes.get();
          var edges = network.body.data.edges.get();
          
          var edgeMap = {};
          edges.forEach(function(e) {
             if (!e.from || !e.to) return;
             var f = String(e.from), t = String(e.to);
             edgeMap[f + '_' + t] = true;
             edgeMap[t + '_' + f] = true;
          });

          var bycat = {};
          nodes.forEach(function(n) {
            if (n.hidden) return;
            var c = String(n[currentAreaAttr] || 'Uncategorized');
            if (c === '' || c === 'NA' || c === 'null' || c === 'undefined') c = 'Uncategorized';
            (bycat[c] || (bycat[c] = [])).push(n.id);
          });
          
          var cats = Object.keys(bycat).filter(function(c) { return c !== 'Uncategorized'; });
          var k = cats.length;
          
          var ncol = Math.max(1, Math.ceil(Math.sqrt(k)));
          var nrow = Math.max(1, Math.ceil(k / ncol));
          var cell = 450; 
          
          var centers = {};
          for (var i = 0; i < k; i++) {
            var r = Math.floor(i / ncol);
            var c = i % ncol;
            centers[cats[i]] = {
              cx: (c - (ncol - 1) / 2) * cell,
              cy: ((nrow - 1) / 2 - r) * cell
            };
          }
          
          var pos = network.getPositions();
          var simNodes = [];
          nodes.forEach(function(n) {
             if (n.hidden) return;
             var cat = String(n[currentAreaAttr] || 'Uncategorized');
             if (cat === '' || cat === 'NA' || cat === 'null' || cat === 'undefined') cat = 'Uncategorized';
             var cx = 0, cy = 0;
             if (centers[cat]) {
                 cx = centers[cat].cx;
                 cy = centers[cat].cy;
             }
             var currentP = pos[n.id] || {x: cx + (Math.random()-0.5)*100, y: cy + (Math.random()-0.5)*100};
             simNodes.push({
                id: String(n.id),
                cat: cat,
                x: currentP.x,
                y: currentP.y,
                vx: 0,
                vy: 0,
                radius: (n.size || 20) + 15 // padding for collision
             });
          });

          // Custom Force Simulation for Organic Clustering
          var iterations = 200;
          var alpha = 1.0;
          for (var iter = 0; iter < iterations; iter++) {
             alpha *= 0.98; // cooling
             
             // 1. Weak attraction to category centers (Uncategorized float freely but weakly pull to center)
             for (var i = 0; i < simNodes.length; i++) {
                var sn = simNodes[i];
                if (sn.cat !== 'Uncategorized' && centers[sn.cat]) {
                   var center = centers[sn.cat];
                   sn.vx += (center.cx - sn.x) * 0.02 * alpha;
                   sn.vy += (center.cy - sn.y) * 0.02 * alpha;
                } else {
                   sn.vx += (0 - sn.x) * 0.003 * alpha;
                   sn.vy += (0 - sn.y) * 0.003 * alpha;
                }
             }
             
             for (var i = 0; i < simNodes.length; i++) {
                for (var j = i + 1; j < simNodes.length; j++) {
                   var n1 = simNodes[i];
                   var n2 = simNodes[j];
                   var dx = n1.x - n2.x;
                   var dy = n1.y - n2.y;
                   var distSq = dx*dx + dy*dy;
                   if (distSq === 0) { dx = Math.random()-0.5; dy = Math.random()-0.5; distSq = dx*dx+dy*dy; }
                   var dist = Math.sqrt(distSq);
                   var minDist = n1.radius + n2.radius;
                   
                   // 2. Edge Attraction
                   if (edgeMap[n1.id + '_' + n2.id]) {
                      var pullEdge = (dist - minDist * 1.5) * 0.008 * alpha;
                      n1.vx -= (dx / dist) * pullEdge;
                      n1.vy -= (dy / dist) * pullEdge;
                      n2.vx += (dx / dist) * pullEdge;
                      n2.vy += (dy / dist) * pullEdge;
                   }
                   
                   // 3. Attraction between nodes of same category (skip Uncategorized)
                   if (n1.cat === n2.cat && n1.cat !== 'Uncategorized') {
                      var pullCat = (dist - minDist * 1.5) * 0.005 * alpha;
                      n1.vx -= (dx / dist) * pullCat;
                      n1.vy -= (dy / dist) * pullCat;
                      n2.vx += (dx / dist) * pullCat;
                      n2.vy += (dy / dist) * pullCat;
                   }
                   
                   // 4. Collision / Repulsion
                   if (dist < minDist) {
                      var force = (minDist - dist) / dist * 0.8 * alpha; // strong collision
                      n1.vx += dx * force;
                      n1.vy += dy * force;
                      n2.vx -= dx * force;
                      n2.vy -= dy * force;
                   } else if (dist < minDist * 4) {
                      var force = (400 / distSq) * alpha; // soft repulsion
                      n1.vx += dx * force;
                      n1.vy += dy * force;
                      n2.vx -= dx * force;
                      n2.vy -= dy * force;
                   }
                }
             }
             
             // 5. Apply velocity and add friction
             for (var i = 0; i < simNodes.length; i++) {
                var sn = simNodes[i];
                sn.x += sn.vx;
                sn.y += sn.vy;
                sn.vx *= 0.6; // strong friction
                sn.vy *= 0.6;
             }
          }

          var updates = simNodes.map(function(sn) {
             return {id: sn.id, x: sn.x, y: sn.y};
          });
          
          network.setOptions({physics: {enabled: false}});
          network.body.data.nodes.update(updates);
          if (typeof refreshNodes === 'function') refreshNodes();
          network.fit({animation: true});
        };
      }

      // Button Container (Flexbox)
      var btnContainer = document.createElement('div');
      Object.assign(btnContainer.style, {
        position: 'absolute', bottom: '15px', right: '15px', zIndex: '1000',
        display: 'flex', flexDirection: 'column', gap: '8px', alignItems: 'center'
      });
      container.appendChild(btnContainer);

      var btnStyle = {
        backgroundColor: 'rgba(255, 255, 255, 0.95)', width: 'clamp(26px, 4vmin, 34px)', height: 'clamp(26px, 4vmin, 34px)',
        borderRadius: '6px', boxShadow: '0 2px 10px rgba(0,0,0,0.1)', border: '1px solid #ddd',
        display: 'flex', alignItems: 'center', justifyContent: 'center', cursor: 'pointer', transition: 'all 0.2s'
      };

      // Settings Button
      var settingsBtn = document.createElement('div');
      Object.assign(settingsBtn.style, btnStyle);
      settingsBtn.innerHTML = \"<svg width='60%' height='60%' viewBox='0 0 24 24' fill='none' stroke='#333' stroke-width='2.5' stroke-linecap='round' stroke-linejoin='round'><circle cx='12' cy='12' r='3'></circle><path d='M19.4 15a1.65 1.65 0 0 0 .33 1.82l.06.06a2 2 0 0 1 0 2.83 2 2 0 0 1-2.83 0l-.06-.06a1.65 1.65 0 0 0-1.82-.33 1.65 1.65 0 0 0-1 1.51V21a2 2 0 0 1-2 2 2 2 0 0 1-2-2v-.09A1.65 1.65 0 0 0 9 19.4a1.65 1.65 0 0 0-1.82.33l-.06.06a2 2 0 0 1-2.83 0 2 2 0 0 1 0-2.83l.06-.06a1.65 1.65 0 0 0 .33-1.82 1.65 1.65 0 0 0-1.51-1H3a2 2 0 0 1-2-2 2 2 0 0 1 2-2h.09A1.65 1.65 0 0 0 4.6 9a1.65 1.65 0 0 0-.33-1.82l-.06-.06a2 2 0 0 1 0-2.83 2 2 0 0 1 2.83 0l.06.06a1.65 1.65 0 0 0 1.82.33H9a1.65 1.65 0 0 0 1-1.51V3a2 2 0 0 1 2-2 2 2 0 0 1 2 2v.09a1.65 1.65 0 0 0 1 1.51 1.65 1.65 0 0 0 1.82-.33l.06-.06a2 2 0 0 1 2.83 0 2 2 0 0 1 0 2.83l-.06.06a1.65 1.65 0 0 0-.33 1.82V9a1.65 1.65 0 0 0 1.51 1H21a2 2 0 0 1 2 2 2 2 0 0 1-2 2h-.09a1.65 1.65 0 0 0-1.51 1z'></path></svg>\";
      settingsBtn.title = \"Ajustes\";
      settingsBtn.onmouseover = function() { this.style.backgroundColor = '#f5f5f5'; };
      settingsBtn.onmouseout = function() { this.style.backgroundColor = 'rgba(255, 255, 255, 0.95)'; };
      settingsBtn.onclick = function() { settingsModal.style.display = (settingsModal.style.display === 'block' ? 'none' : 'block'); };
      // Info Modal and Button are defined below, we'll append Settings after Info.
      // Info Modal
      var infoModal = document.createElement('div');
      infoModal.id = 'digraph_info_modal';
      Object.assign(infoModal.style, {
        position: 'absolute', top: '50%', left: '50%', transform: 'translate(-50%, -50%)',
        width: '90%', maxWidth: '450px', maxHeight: '80vh', overflowY: 'auto', boxSizing: 'border-box',
        backgroundColor: '#fff', zIndex: '2000', padding: '20px', borderRadius: '8px', 
        boxShadow: '0 4px 20px rgba(0,0,0,0.2)', border: '1px solid #eaeaea', display: 'none', 
        fontFamily: 'Inter, Roboto, sans-serif'
      });
      infoModal.innerHTML = '<div style=\"display:flex; justify-content:space-between; align-items:center; border-bottom:1px solid #eaeaea; padding-bottom:10px; margin-bottom:15px;\">' +
                            '<h3 style=\"margin:0; color:#444; font-size:16px;\">' + (x.sim_data ? (x.dict.network_view || 'Vista de Red') : x.dict.self_digraph) + '</h3>' +
                            '<span id=\"close_info_modal\" style=\"cursor:pointer; font-size:20px; font-weight:bold; color:#888; line-height:1;\">&times;</span>' +
                            '</div>' +
                            '<p style=\"margin:0; color:#666; font-size:13px; line-height:1.6;\">' + (x.sim_data ? (x.dict.info_text_sim_network || 'Simulación.') : x.dict.info_text_digraph) + '</p>';
      container.appendChild(infoModal);

      // Info Button
      var infoPanel = document.createElement('div');
      Object.assign(infoPanel.style, btnStyle);
      infoPanel.innerHTML = \"<svg width='60%' height='60%' viewBox='0 0 24 24' fill='none' stroke='#333' stroke-width='2.5' stroke-linecap='round' stroke-linejoin='round'><circle cx='12' cy='12' r='10'></circle><line x1='12' y1='16' x2='12' y2='12'></line><line x1='12' y1='8' x2='12.01' y2='8'></line></svg>\";
      infoPanel.title = x.dict.info;
      infoPanel.onmouseover = function() { this.style.backgroundColor = '#f5f5f5'; };
      infoPanel.onmouseout = function() { this.style.backgroundColor = 'rgba(255, 255, 255, 0.95)'; };
      infoPanel.onclick = function() { infoModal.style.display = (infoModal.style.display === 'block' ? 'none' : 'block'); };
      btnContainer.appendChild(infoPanel);
      btnContainer.appendChild(settingsBtn);
      
      infoModal.querySelector('#close_info_modal').onclick = function() { infoModal.style.display = 'none'; };

      // Minimalist Export Panel
      var exportPanel = document.createElement('div');
      Object.assign(exportPanel.style, btnStyle);
      exportPanel.innerHTML = \"<svg width='60%' height='60%' viewBox='0 0 24 24' fill='none' stroke='#333' stroke-width='2.5' stroke-linecap='round' stroke-linejoin='round'><path d='M21 15v4a2 2 0 0 1-2 2H5a2 2 0 0 1-2-2v-4'></path><polyline points='7 10 12 15 17 10'></polyline><line x1='12' y1='15' x2='12' y2='3'></line></svg>\";
      exportPanel.title = x.dict.export_png;
      exportPanel.onmouseover = function() { this.style.backgroundColor = '#f5f5f5'; };
      exportPanel.onmouseout = function() { this.style.backgroundColor = 'rgba(255, 255, 255, 0.95)'; };
      exportPanel.onclick = exportPNG;
      btnContainer.appendChild(exportPanel);

      // Fullscreen Button
      var fsPanel = document.createElement('div');
      Object.assign(fsPanel.style, btnStyle);
      fsPanel.innerHTML = \"<svg width='60%' height='60%' viewBox='0 0 24 24' fill='none' stroke='#333' stroke-width='2.5' stroke-linecap='round' stroke-linejoin='round'><path d='M8 3H5a2 2 0 0 0-2 2v3m18 0V5a2 2 0 0 0-2-2h-3m0 18h3a2 2 0 0 0 2-2v-3M3 16v3a2 2 0 0 0 2 2h3'></path></svg>\";
      fsPanel.title = x.dict.fullscreen;
      fsPanel.onmouseover = function() { this.style.backgroundColor = '#f5f5f5'; };
      fsPanel.onmouseout = function() { this.style.backgroundColor = 'rgba(255, 255, 255, 0.95)'; };
      fsPanel.onclick = function() {
        var fullEl = container.closest('.wt-tab-container') || container;
        if (!document.fullscreenElement) {
          fullEl.requestFullscreen().catch(err => console.log('Error full screen:', err));
        } else {
          document.exitFullscreen();
        }
      };
      btnContainer.appendChild(fsPanel);

      // PCSD Chart Overlay
      var chartContainer = document.createElement('div');
      chartContainer.id = 'pcsd_chart_container';
      Object.assign(chartContainer.style, {
        position: 'absolute', top: '0', left: '0', width: '100%', height: '100%',
        backgroundColor: '#fff', zIndex: '5', display: 'none', padding: '60px 40px 110px 40px', boxSizing: 'border-box'
      });
      chartContainer.innerHTML = '<canvas id=\"pcsd_canvas\"></canvas>';
      container.appendChild(chartContainer);

      // --- Visualization Content ---
      var visHTML = '<div style=\"margin-bottom:15px;\">' +
                    '<label style=\"display:block; margin-bottom:5px; font-weight:bold; color:#444;\">' + x.dict.color_palette + '</label>' +
                    '<select id=\"palette_sel\" style=\"width:100%; padding:4px; border-radius:4px;\">' +
                    Object.keys(x.color_palette_js).map(k => '<option value=\"' + k + '\"' + (k === x.initial_palette ? ' selected' : '') + '>' + ((x.dict.palette_labels && x.dict.palette_labels[k]) || k.charAt(0).toUpperCase() + k.slice(1)) + '</option>').join('') +
                    '</select></div>';
      
      visHTML += '<div style=\"margin-bottom:15px;\">' +
                 '<label style=\"display:block; margin-bottom:5px; font-weight:bold; color:#444;\">' + x.dict.layout_algo + '</label>' +
                 '<select id=\"layout_sel\" style=\"width:100%; padding:4px; border-radius:4px;\">' +
                 Object.keys(x.layouts).map(k => '<option value=\"' + k + '\"' + (k === x.initial_layout ? ' selected' : '') + '>' + ((x.dict.layout_labels && x.dict.layout_labels[k]) || k.charAt(0).toUpperCase() + k.slice(1)) + '</option>').join('') +
                 '</select></div>';
      
      visHTML += '<div style=\"margin-bottom:15px; border-top:1px solid #eee; padding-top:10px;\">' +
                 '<div style=\"display:flex; justify-content:space-between; margin-bottom:5px;\"><b style=\"color:#444; font-size:13px;\">' + x.dict.edge_opacity + '</b></div>' +
                 '<input type=\"range\" id=\"opacity_slider\" min=\"0.1\" max=\"1\" step=\"0.1\" value=\"0.5\" style=\"width:100%; accent-color:#8cc63f;\">' +
                 
                 '<div style=\"display:flex; justify-content:space-between; margin-bottom:5px; margin-top:12px;\"><b style=\"color:#444; font-size:13px;\">' + x.dict.edge_filter + '</b><span id=\"weight_val_txt\" style=\"font-family:monospace; font-size:12px;\">0.00</span></div>' +
                 '<input type=\"range\" id=\"weight_slider\" min=\"0\" max=\"' + x.max_weight + '\" step=\"0.01\" value=\"0\" style=\"width:100%;\">' +
                 '<label style=\"display:flex; align-items:center; margin-top:8px; font-size:11px; cursor:pointer; color:#555;\">' +
                 '<input type=\"checkbox\" id=\"hide_direct_check\" ' + (x.hide_direct ? 'checked' : '') + ' style=\"margin-right:6px;\"> ' + x.dict.hide_direct + '</label>' +
                 '</div>';

      visHTML += '<div style=\"margin-bottom:15px;\">' +
                 '<div style=\"display:flex; justify-content:space-between; margin-bottom:5px;\"><b style=\"color:#444; font-size:13px;\">' + x.dict.node_size + '</b></div>' +
                 '<input type=\"range\" id=\"size_slider\" min=\"0.5\" max=\"3\" step=\"0.1\" value=\"1\" style=\"width:100%; accent-color:#8cc63f;\">' +
                 '</div>';

      visHTML += '<div style=\"margin-bottom:15px;\">' +
                 '<div style=\"display:flex; justify-content:space-between; margin-bottom:5px;\"><b style=\"color:#444; font-size:13px;\">' + x.dict.text_size + '</b></div>' +
                 '<input type=\"range\" id=\"text_size_slider\" min=\"10\" max=\"40\" step=\"1\" value=\"20\" style=\"width:100%; accent-color:#8cc63f;\">' +
                 '</div>';

      visHTML += '<div style=\"margin-bottom:15px; border-top:1px solid #f0f0f0; padding-top:10px; display:flex; gap:5px;\">' +
                 '<button id=\"eraser_tool\" title=\"Click nodes/edges to hide them\" style=\"flex:1; padding:6px; background:#fff; border:1px solid #ccc; border-radius:4px; cursor:pointer; font-size:11px; display:flex; align-items:center; justify-content:center; gap:4px;\">🧹 ' + x.dict.eraser_mode + '</button>' +
                 '<button id=\"btn_reset\" title=\"Reset all settings and restore elements\" style=\"flex:1; padding:6px; background:#f8f9fa; border:1px solid #ccc; border-radius:4px; font-weight:bold; color:#555; cursor:pointer; font-size:11px; display:flex; align-items:center; justify-content:center; gap:4px;\">↺ ' + x.dict.reset + '</button>' +
                 '</div>';

      visHTML += '<div style=\"border-top:1px solid #eee; padding-top:10px;\">' +
                 '<b style=\"color:#444; display:block; margin-bottom:8px; font-size:13px;\">' + x.dict.visible_constructs + '</b>' +
                 '<div style=\"display:flex; gap:5px; margin-bottom:8px;\">' +
                 '<button id=\"sel_all\" style=\"flex:1; font-size:10px; cursor:pointer;\">' + x.dict.all + '</button>' +
                 '<button id=\"sel_none\" style=\"flex:1; font-size:10px; cursor:pointer;\">' + x.dict.none_btn + '</button></div>' +
                 '<div id=\"node_list\" style=\"max-height:150px; overflow-y:auto; border:1px solid #f0f0f0; padding:5px;\"></div></div>';
      
      visContent.innerHTML = visHTML;
      settingsModalBody.appendChild(visContent);

      // --- Simulation Panel ---
      if (x.sim_data) {
        var sim = x.sim_data;
        network._simHistory = [];
        network._simCurrentI = 0;
        var targetSelf = [...sim.initial_self];
        var simMaxIter = sim.max_iter || 10;
        var simThreshold = sim.threshold || 'saturation';

        // ── 1. Settings Panel – collapsible, Bottom-Left ──────────────────
        var simSettingsContent = createPanel('sim_settings_panel', x.dict.sim_settings,
          {bottom: '10px', left: '10px', width: '260px', borderLeft: '4px solid #3498db'});
        
        var settingsHTML =
          (isSidebarMode ? '' :
          '<div style=\"margin-bottom:12px; border-bottom:1px solid #eee; padding-bottom:10px; display:flex; gap:5px;\">' +
            '<button id=\"btn_view_graph\" style=\"flex:1; padding:6px; background:#f4f9ef; border:1px solid #8cc63f; border-radius:4px; font-weight:bold; color:#5c8822; cursor:pointer; font-size:11px;\">' + x.dict.network_view + '</button>' +
            '<button id=\"btn_view_pcsd\" style=\"flex:1; padding:6px; background:#fff; border:1px solid #ccc; border-radius:4px; font-weight:bold; color:#555; cursor:pointer; font-size:11px;\">' + x.dict.pcsd_chart + '</button>' +
          '</div>') +
          '<div style=\"margin-bottom:12px;\">' +
            '<label style=\"display:block; margin-bottom:5px; font-weight:bold; color:#444;\">' + x.dict.thr_function + '</label>' +
            '<select id=\"sim_thr_sel\" style=\"width:100%; padding:4px; border-radius:4px;\">' +
              '<option value=\"saturation\"' + (simThreshold === 'saturation' ? ' selected' : '') + '>' + x.dict.saturation + '</option>' +
              '<option value=\"tanh\"' + (simThreshold === 'tanh' ? ' selected' : '') + '>' + x.dict.tanh + '</option>' +
              '<option value=\"linear\"' + (simThreshold === 'linear' ? ' selected' : '') + '>' + x.dict.linear + '</option>' +
            '</select>' +
          '</div>' +
          '<div style=\"margin-bottom:12px;\">' +
            '<label style=\"display:block; margin-bottom:5px; font-weight:bold; color:#444;\">' + x.dict.sim_depth + '</label>' +
            '<input type=\"number\" id=\"sim_depth_input\" min=\"1\" max=\"50\" value=\"' + simMaxIter + '\"' +
              ' style=\"width:100%; padding:4px; border-radius:4px; border:1px solid #ccc;\">' +
          '</div>' +
          '<div style=\"margin-bottom:12px; border-top:1px solid #eee; padding-top:10px;\">' +
            '<label style=\"display:block; margin-bottom:5px; font-weight:bold; color:#444;\">' + x.dict.playback_speed + '</label>' +
            '<input type=\"range\" id=\"speed_slider\" min=\"0.25\" max=\"2\" step=\"0.25\" value=\"1\" style=\"width:100%; accent-color:#8cc63f;\">' +
            '<div id=\"speed_txt\" style=\"text-align:right; font-size:10px; color:#888; margin-top:2px;\">1.00x</div>' +
          '</div>' +
          '<div style=\"border-top:1px solid #eee; padding-top:10px; display:flex; flex-direction:column; flex:1; min-height:0;\">' +
            '<b style=\"color:#444; display:block; margin-bottom:8px; flex-shrink:0;\">' + x.dict.scenario + '</b>' +
            '<div id=\"act_list\" style=\"flex:1; overflow-y:auto; overflow-x:hidden; border:1px solid #f8f9fa; padding:5px; background:#fafafa;\"></div>' +
          '</div>' +
          '<div style=\"border-top:1px solid #eee; padding-top:10px; margin-top:8px; flex-shrink:0;\">' +
            '<button id=\"reset_sim\" style=\"width:100%; padding:5px 0; background:#fef9f1; border:1px solid #e0c97a; border-radius:4px; font-size:11px; font-weight:bold; color:#888; cursor:pointer;\">&#8635; ' + x.dict.reset_scenario + '</button>' +
          '</div>';

        simSettingsContent.innerHTML = settingsHTML;

        // ── 2. Timeline – minimal single-line pill ────────────────────────
        var timelineEl = document.createElement('div');
        timelineEl.id = 'timeline_panel';
        Object.assign(timelineEl.style, {
          position: 'absolute', zIndex: '1000',
          backgroundColor: 'rgba(255,255,255,0.95)',
          padding: '7px 10px', borderRadius: '20px',
          boxShadow: '0 2px 10px rgba(0,0,0,0.12)',
          border: '1px solid #d0e8f8',
          fontFamily: 'Segoe UI, Tahoma, sans-serif', fontSize: '12px',
          bottom: '20px', left: '50%', transform: 'translateX(-50%)', width: '330px', boxSizing: 'border-box'
        });
        timelineEl.innerHTML =
          '<div style=\"display:flex; align-items:center; gap:7px; width:100%;\">' +
            '<button id=\"play_btn\" title=\"Play\" style=\"min-width:26px;min-height:26px;width:26px;height:26px;flex-shrink:0;border:none;border-radius:50%;background:#8cc63f;color:#fff;font-size:12px;cursor:pointer;padding:0;display:flex;align-items:center;justify-content:center;box-sizing:border-box;\">&#9654;</button>' +
            '<button id=\"pause_btn\" title=\"Pause\" style=\"min-width:26px;min-height:26px;width:26px;height:26px;flex-shrink:0;border:1px solid #ccc;border-radius:50%;background:#f5f5f5;color:#555;font-size:10px;cursor:pointer;padding:0;display:flex;align-items:center;justify-content:center;box-sizing:border-box;\">&#9646;&#9646;</button>' +
            '<button id=\"stop_btn\" title=\"Stop (Reset)\" style=\"min-width:26px;min-height:26px;width:26px;height:26px;flex-shrink:0;border:1px solid #ccc;border-radius:50%;background:#f5f5f5;color:#e74c3c;font-size:12px;cursor:pointer;padding:0;display:flex;align-items:center;justify-content:center;box-sizing:border-box;\">&#9632;</button>' +
            '<input type=\"range\" id=\"sim_slider\" min=\"0\" max=\"' + simMaxIter + '\" value=\"0\"' +
              ' style=\"flex:1;accent-color:#8cc63f;cursor:pointer;margin:0;\">' +
            '<span id=\"iter_label\" style=\"flex-shrink:0;font-size:11px;font-weight:bold;color:#8cc63f;white-space:nowrap;min-width:38px;text-align:right;\">0/' + simMaxIter + '</span>' +
          '</div>';
        container.appendChild(timelineEl);

        // ── Activation sliders ────────────────────────────────────────────
        var actList = simSettingsContent.querySelector('#act_list');

        // Inject slider reset CSS once (remove native left-fill)
        (function() {
          if(document.getElementById('wt_slider_css')) return;
          var s = document.createElement('style');
          s.id = 'wt_slider_css';
          s.textContent =
            '.wt-slider{-webkit-appearance:none;appearance:none;width:100%;height:6px;background:transparent;outline:none;cursor:pointer;padding:0;margin:4px 0;}' +
            '.wt-slider::-webkit-slider-thumb{-webkit-appearance:none;width:13px;height:13px;border-radius:50%;background:#555;cursor:pointer;margin-top:-4px;border:2px solid #fff;box-shadow:0 0 2px rgba(0,0,0,0.3);}' +
            '.wt-slider::-moz-range-thumb{width:13px;height:13px;border-radius:50%;background:#555;border:2px solid #fff;cursor:pointer;}' +
            '.wt-slider::-webkit-slider-runnable-track{height:6px;background:var(--wt-track-bg,#e0e0e0);border-radius:3px;}' +
            '.wt-slider::-moz-range-track{height:6px;background:var(--wt-track-bg,#e0e0e0);border-radius:3px;}' +
            'input[type=\"range\"], input[type=\"checkbox\"], input[type=\"radio\"] { accent-color: #8cc63f !important; }' +
            'select:focus { border-color: #8cc63f !important; outline: none; }' +
            'option:checked { background-color: #8cc63f !important; color: white !important; }' +
            'option:hover, option:focus, option:active { background-color: #8cc63f !important; color: white !important; box-shadow: 0 0 10px 100px #8cc63f inset !important; }' +
            '::selection { background-color: #8cc63f !important; color: white !important; }';
          document.head.appendChild(s);
        })();

        // Helper: compute fill gradient string
        var getSliderFill = function(initVal, curVal, idealVal) {
          var hasIdeal = (idealVal !== 0);
          var delta = curVal - initVal;
          var color;
          if(!hasIdeal) {
            color = (Math.abs(delta) < 0.001) ? null : '#F0C040';
          } else {
            var towardIdeal = Math.sign(curVal - initVal) === Math.sign(idealVal - initVal);
            color = (Math.abs(delta) < 0.001) ? null : (towardIdeal ? '#4CAF50' : '#E53935');
          }
          if(!color) return '#e0e0e0';
          var initPct = ((initVal + 1) / 2) * 100;
          var curPct  = ((curVal  + 1) / 2) * 100;
          var left    = Math.min(initPct, curPct).toFixed(1);
          var right   = Math.max(initPct, curPct).toFixed(1);
          return 'linear-gradient(to right,#e0e0e0 ' + left + '%,' + color + ' ' + left + '%,' + color + ' ' + right + '%,#e0e0e0 ' + right + '%)';
        };

        // Set gradient via CSS custom property (pseudo-elements inherit it)
        var updateSliderTrack = function(sliderEl, initVal, idealVal) {
          sliderEl.style.setProperty('--wt-track-bg', getSliderFill(initVal, parseFloat(sliderEl.value), idealVal));
        };

        sim.lpoles.forEach(function(lp, i) {
          var initVal  = sim.initial_self[i];
          var idealVal = (sim.initial_ideal && sim.initial_ideal[i] !== undefined) ? sim.initial_ideal[i] : 0;

          var div = document.createElement('div');
          div.style.marginBottom = '12px';
          div.innerHTML =
            '<div style=\"display:flex; justify-content:space-between; font-size:10px; color:#555; margin-bottom:3px;\">' +
              '<span style=\"white-space:nowrap; overflow:hidden; text-overflow:ellipsis; width:150px;\" title=\"' + lp + ' - ' + sim.rpoles[i] + '\">' + lp + ' - ' + sim.rpoles[i] + '</span>' +
              '<b class=\"val-badge\" style=\"color:#555;\">' + initVal.toFixed(2) + '</b>' +
            '</div>' +
            '<input type=\"range\" class=\"target-slider wt-slider\" data-idx=\"' + i + '\"' +
              ' data-init=\"' + initVal + '\" data-ideal=\"' + idealVal + '\"' +
              ' min=\"-1\" max=\"1\" step=\"0.05\" value=\"' + initVal + '\">';
          actList.appendChild(div);

          // Apply initial track style
          var sliderEl = div.querySelector('.wt-slider');
          updateSliderTrack(sliderEl, initVal, idealVal);
        });


        // ── Simulation engine ─────────────────────────────────────────────
        var runSimulation = function() {
          var m = parseInt(simSettingsContent.querySelector('#sim_depth_input').value) || simMaxIter;
          var thr_type = simSettingsContent.querySelector('#sim_thr_sel').value;
          var thr = function(v) {
            return thr_type === 'tanh' ? Math.tanh(v) :
                   thr_type === 'saturation' ? Math.max(-1, Math.min(1, v)) : v;
          };
          var s = [...sim.initial_self];
          var a = targetSelf.map(function(t, idx) { return t - sim.initial_self[idx]; });
          var h = [[...s]];
          for(var i = 0; i < m; i++) {
            var next = s.map(function(v, idx) { return thr(v + a[idx]); });
            var delta = next.map(function(v, idx) { return v - s[idx]; });
            s = next;
            var newA = new Array(s.length).fill(0);
            for(var r = 0; r < s.length; r++)
              for(var c = 0; c < s.length; c++)
                newA[r] += sim.weights[c][r] * delta[c];
            a = newA;
            h.push([...s]);
          }
          network._simHistory = h;
          // Sync slider range
          var sliderEl = timelineEl.querySelector('#sim_slider');
          sliderEl.max = m;
          if(network._simCurrentI > m) network._simCurrentI = m;
          refreshNodes();
          if(typeof updatePcsdData === 'function') updatePcsdData();
        };

        var updateIteration = function(idx) {
          network._simCurrentI = idx;
          if(!network._simHistory[idx]) return;
          timelineEl.querySelector('#sim_slider').value = idx;
          timelineEl.querySelector('#iter_label').innerText = idx + '/' + (parseInt(simSettingsContent.querySelector('#sim_depth_input').value) || simMaxIter);
          refreshNodes();
        };

        // ── Event listeners ───────────────────────────────────────────────
        simSettingsContent.querySelector('#sim_thr_sel').onchange = runSimulation;
        simSettingsContent.querySelector('#sim_depth_input').onchange = runSimulation;

        var playbackSpeed = 1.0;
        simSettingsContent.querySelector('#speed_slider').oninput = function() {
          playbackSpeed = parseFloat(this.value);
          if (playbackSpeed <= 0) playbackSpeed = 0.25;
          simSettingsContent.querySelector('#speed_txt').innerText = playbackSpeed.toFixed(2) + 'x';
        };

        timelineEl.querySelector('#sim_slider').oninput = function() { updateIteration(parseInt(this.value)); };

        var simTimer = null;
        var _playStep = function(m) {
          if(network._simCurrentI >= m) { simTimer = null; return; }
          var from = network._simCurrentI;
          var to   = from + 1;
          
          timelineEl.querySelector('#sim_slider').value = to;
          timelineEl.querySelector('#iter_label').innerText = to + '/' + m;
          
          var tweenDur = 1500 / playbackSpeed;
          var stepInt  = 1800 / playbackSpeed;
          
          tweenToIteration(from, to, tweenDur);
          network._simCurrentI = to;
          simTimer = setTimeout(function() { _playStep(m); }, stepInt);
        };

        timelineEl.querySelector('#play_btn').onclick = function() {
          if(simTimer) return;
          var m = parseInt(simSettingsContent.querySelector('#sim_depth_input').value) || simMaxIter;
          if(network._simCurrentI >= m) { network._simCurrentI = 0; refreshNodes(); }
          _playStep(m);
        };
        timelineEl.querySelector('#pause_btn').onclick = function() {
          if(simTimer) { clearTimeout(simTimer); simTimer = null; }
          if(_tweenRAF) { cancelAnimationFrame(_tweenRAF); _tweenRAF = null; }
        };
        timelineEl.querySelector('#stop_btn').onclick = function() {
          if(simTimer) { clearTimeout(simTimer); simTimer = null; }
          if(_tweenRAF) { cancelAnimationFrame(_tweenRAF); _tweenRAF = null; }
          network._activePulses = [];
          updateIteration(0);
        };
        simSettingsContent.querySelector('#reset_sim').onclick = function() {
          if(simTimer) clearTimeout(simTimer);
          targetSelf = [...sim.initial_self];
          actList.querySelectorAll('.target-slider').forEach(function(s, idx) {
            var initV  = parseFloat(s.getAttribute('data-init'));
            var idealV = parseFloat(s.getAttribute('data-ideal'));
            s.value = targetSelf[idx];
            s.parentNode.querySelector('.val-badge').innerText = targetSelf[idx].toFixed(2);
            updateSliderTrack(s, initV, idealV);
          });
          runSimulation();
          updateIteration(0);
          if(typeof updatePcsdData === 'function') updatePcsdData();
        };
        actList.oninput = function(e) {
          if(e.target.classList.contains('target-slider')) {
            var idx = parseInt(e.target.getAttribute('data-idx'));
            var val = parseFloat(e.target.value);
            var initV  = parseFloat(e.target.getAttribute('data-init'));
            var idealV = parseFloat(e.target.getAttribute('data-ideal'));
            targetSelf[idx] = val;
            e.target.parentNode.querySelector('.val-badge').innerText = val.toFixed(2);
            updateSliderTrack(e.target, initV, idealV);
            runSimulation();
          }
        };

        // --- PCSD View & Chart Logic ---
        var btnGraph = simSettingsContent.querySelector('#btn_view_graph');
        var btnPcsd = simSettingsContent.querySelector('#btn_view_pcsd');
        var plotlyCanvasId = isSidebarMode ? 'wsim_pcsd_plot' : 'pcsd_canvas_plotly';
        var _pcsdInit = false;

        if (btnGraph && btnPcsd) {
          var updateView = function(mode) {
            var isGraph = (mode === 'graph');
            chartContainer.style.display = isGraph ? 'none' : 'block';
            btnGraph.style.background = isGraph ? '#f4f9ef' : '#fff';
            btnGraph.style.border = isGraph ? '1px solid #8cc63f' : '1px solid #ccc';
            btnGraph.style.color = isGraph ? '#5c8822' : '#555';
            
            btnPcsd.style.background = !isGraph ? '#f4f9ef' : '#fff';
            btnPcsd.style.border = !isGraph ? '1px solid #8cc63f' : '1px solid #ccc';
            btnPcsd.style.color = !isGraph ? '#5c8822' : '#555';
            if (!isGraph) {
              if (window.Plotly && !_pcsdInit) {
                initPcsdChart();
              }
            }
          };
          btnGraph.onclick = function() { updateView('graph'); };
          btnPcsd.onclick = function() { updateView('pcsd'); };
        } else {
           // If we are in sidebar mode, just initialize it silently so Plotly reacts
           setTimeout(function() { if(window.Plotly && !_pcsdInit) initPcsdChart(); }, 500);
        }

        var _plotlyLayout = {
          margin: { l: 60, r: 20, t: 30, b: 50 },
          showlegend: true,
          legend: { title: { text: '<b>PERSONAL CONSTRUCTS</b>' }, font: {size: 11} },
          xaxis: { title: 'ITERATIONS', zeroline: false },
          yaxis: { title: 'SELF DIFFERENTIAL', zeroline: true, zerolinecolor: '#666', zerolinewidth: 2, range: [-2, 2] }
        };

        var initPcsdChart = function() {
          if (!window.Plotly) return;
          var poles = sim.lpoles.map((lp, i) => lp + ' - ' + sim.rpoles[i]);
          var plotlyPalette = [
            '#636EFA', '#EF553B', '#00CC96', '#AB63FA', '#FFA15A', '#19D3F3', '#FF6692', '#B6E880', 
            '#FF97FF', '#FECB52', '#0d0887', '#46039f', '#7201a8', '#9c179e', '#bd3786', '#d8576b'
          ];
          var symbols = ['circle', 'square', 'diamond', 'cross', 'x', 'triangle-up'];
          
          var traces = poles.map((p, i) => {
            return {
              name: (p.length > 30 ? p.substring(0, 27) + '...' : p),
              x: [],
              y: [],
              mode: 'lines+markers',
              line: { color: plotlyPalette[i % plotlyPalette.length], width: 2.5, shape: 'spline' },
              marker: { symbol: symbols[Math.floor(i / plotlyPalette.length) % symbols.length], size: 8 }
            };
          });

          var canvasEl = document.getElementById(plotlyCanvasId);
          if (!canvasEl) {
              console.warn('Plotly canvas element not found: ' + plotlyCanvasId);
              return;
          }
          _pcsdInit = true;
          // Clear any R htmlwidget content first to avoid conflicts
          canvasEl.innerHTML = '';
          try {
              Plotly.newPlot(plotlyCanvasId, traces, _plotlyLayout, {responsive: true, displayModeBar: false});
          } catch (err) {
              console.error('Plotly.newPlot error:', err);
          }
          updatePcsdData();
        };

        var updatePcsdData = function() {
          if (!window.Plotly || !_pcsdInit || !network._simHistory) return;
          var history = network._simHistory;
          var initial = sim.initial_self;
          
          var x_data = history.map((_, i) => i);
          
          var tracesUpdate = {
            x: Array(initial.length).fill(x_data),
            y: initial.map((_, dsIdx) => history.map(step => step[dsIdx] - initial[dsIdx]))
          };
          
          Plotly.update(plotlyCanvasId, tracesUpdate, {});
        };

        runSimulation();
      }

      visContent.querySelector('#btn_reset').onclick = function() {
        visContent.querySelector('#opacity_slider').value = 0.5;
        visContent.querySelector('#size_slider').value = 1.0;
        visContent.querySelector('#text_size_slider').value = 20;
        currentDistMult = 1.0;
        currentSizeMult = 1.0;
        currentTextSize = 20;
        
        // Restore Visibility
        var nodesDS = network.body.data.nodes;
        var edgesDS = network.body.data.edges;
        nodesDS.update(nodesDS.get().map(n => ({id: n.id, hidden: false})));
        edgesDS.update(edgesDS.get().map(e => ({id: e.id, hidden: false})));
        visContent.querySelectorAll('.node-check').forEach(chk => chk.checked = true);
        

        network.unselectAll();

        refreshNodes();
        network.stabilize();
        if(typeof initPcsdChart === 'function') initPcsdChart();
      };

      visContent.querySelector('#palette_sel').onchange = function() { refreshNodes(); };
      visContent.querySelector('#layout_sel').onchange = function() {
        var layout = x.layouts[this.value];
        if(!layout) return;
        network.setOptions({physics: {enabled: false}});
        var updates = layout.id.map((id, i) => ({id: String(id), x: layout.x[i] * currentDistMult, y: layout.y[i] * currentDistMult}));
        network.body.data.nodes.update(updates);
        refreshNodes();
        network.fit({animation: true});
      };

      visContent.querySelector('#opacity_slider').oninput = function() {
        refreshNodes();
      };

      visContent.querySelector('#size_slider').oninput = function() {
        currentSizeMult = parseFloat(this.value);
        refreshNodes();
      };

      visContent.querySelector('#text_size_slider').oninput = function() {
        currentTextSize = parseInt(this.value);
        refreshNodes();
      };

      var updateEdges = function() {
        var threshold = parseFloat(visContent.querySelector('#weight_slider').value);
        var hideDirect = visContent.querySelector('#hide_direct_check').checked;
        visContent.querySelector('#weight_val_txt').innerText = threshold.toFixed(2);
        var edgesDS = network.body.data.edges;
        var updates = edgesDS.get().map(edge => ({
          id: edge.id, 
          hidden: Math.abs(edge.weight) < threshold || (hideDirect && (edge.weight > 0 || edge.is_dilemmatic))
        }));
        edgesDS.update(updates);
      };

      visContent.querySelector('#weight_slider').oninput = updateEdges;
      visContent.querySelector('#hide_direct_check').onchange = updateEdges;

      // ── Custom Network Tools Logic ─────────────────────────────────────────
      var eraserActive = false;
      var eraserBtn = visContent.querySelector('#eraser_tool');
      
      eraserBtn.onclick = function() {
        eraserActive = !eraserActive;
        this.style.background = eraserActive ? '#fff0f0' : '#fff';
        this.style.borderColor = eraserActive ? '#e74c3c' : '#ccc';
        this.style.color = eraserActive ? '#e74c3c' : '#333';
        network.canvas.body.container.style.cursor = eraserActive ? 'crosshair' : 'default';
      };




      var nodeClickState = {};
      
      function hex2rgba(hex, alpha) {
          if(hex && hex.startsWith('#')) {
              var r = parseInt(hex.slice(1, 3), 16),
                  g = parseInt(hex.slice(3, 5), 16),
                  b = parseInt(hex.slice(5, 7), 16);
              return 'rgba(' + (isNaN(r)?150:r) + ',' + (isNaN(g)?150:g) + ',' + (isNaN(b)?150:b) + ',' + alpha + ')';
          }
          return hex; 
      }

      network.on('click', function(params) {
        if(eraserActive) {
          if(params.nodes.length > 0) {
            var nodeId = params.nodes[0];
            network.body.data.nodes.update({id: nodeId, hidden: true});
            var chk = visContent.querySelector('input[value=\"' + nodeId + '\"]');
            if(chk) chk.checked = false;
          } else if(params.edges.length > 0) {
            network.body.data.edges.update({id: params.edges[0], hidden: true});
          }
          return;
        }
        
        var edgesDS = network.body.data.edges;
        var nodesDS = network.body.data.nodes;
        var scheme = visContent.querySelector('#palette_sel').value;

        if (params.nodes.length === 1) {
          var nodeId = params.nodes[0];
          var state = nodeClickState[nodeId] || 0;
          state = (state + 1) % 4; // 1: OUT, 2: IN, 3: BOTH, 0: NONE -> skip 0
          if (state === 0) state = 1;
          
          nodeClickState = {};
          nodeClickState[nodeId] = state;

          var connectedNodes = {};
          connectedNodes[nodeId] = true;

          var edgeUpdates = edgesDS.get().map(function(e) {
            var isOut = String(e.from) === String(nodeId);
            var isIn = String(e.to) === String(nodeId);
            
            var highlight = false;
            if (state === 1 && isOut) highlight = true;
            if (state === 2 && isIn) highlight = true;
            if (state === 3 && (isIn || isOut)) highlight = true;

            if (highlight) {
               connectedNodes[e.from] = true;
               connectedNodes[e.to] = true;
            }

            var baseColor = e.orig_color;
            var highColor = e.orig_highlight;
            if (scheme === 'grey scale') {
               baseColor = 'rgba(153,153,153,0.5)';
               highColor = 'rgba(153,153,153,1)';
            }

            if (highlight) {
              return {id: e.id, color: {color: highColor, highlight: highColor}, hidden: false};
            } else {
              return {id: e.id, color: {color: 'rgba(200,200,200,0.05)', highlight: 'rgba(200,200,200,0.05)'}};
            }
          });
          
          refreshNodes();
          edgesDS.update(edgeUpdates);

          var nodeUpdates = nodesDS.get().map(function(n) {
             if (connectedNodes[n.id]) {
                if (String(n.id) === String(nodeId)) {
                   // Ensure the selected node remains exactly its native color without any highlight darkening
                   return {
                      id: n.id,
                      color: {
                         background: n.color.background,
                         border: n.color.border,
                         highlight: { background: n.color.background, border: n.color.border }
                      },
                      font: { color: '#000000', strokeColor: '#ffffff' }
                   };
                } else {
                   return {
                      id: n.id,
                      color: {
                         background: hex2rgba(n.color.background, 0.4),
                         border: hex2rgba(n.color.border, 0.4),
                         highlight: { background: hex2rgba(n.color.background, 0.4), border: hex2rgba(n.color.border, 0.4) }
                      },
                      font: { color: '#000000', strokeColor: '#ffffff' }
                   };
                }
             } else {
                return {
                   id: n.id,
                   color: {
                      background: 'rgba(200,200,200,0.1)',
                      border: 'rgba(200,200,200,0.1)',
                      highlight: { background: 'rgba(200,200,200,0.1)', border: 'rgba(200,200,200,0.1)' }
                   },
                   font: { color: 'rgba(0,0,0,0)', strokeColor: 'rgba(0,0,0,0)' }
                };
             }
          });
          nodesDS.update(nodeUpdates);
          
        } else if (params.nodes.length === 0) {
          nodeClickState = {};
          refreshNodes();
        }
      });

      // Initial refresh to ensure centering on load
      refreshNodes();

      var nodeListDiv = visContent.querySelector('#node_list');
      var nodesDS = network.body.data.nodes;
      nodesDS.get().forEach(node => {
        var div = document.createElement('div');
        div.className = 'node-item';
        div.style.marginBottom = '4px';
        div.style.display = 'flex';
        div.style.alignItems = 'center';
        div.innerHTML = '<label style=\"cursor:pointer; display:flex; align-items:center; font-size:11px;\">' +
                        '<input type=\"checkbox\" class=\"node-check\" value=\"' + node.id + '\" data-id=\"' + node.id + '\" ' + (node.hidden ? '' : 'checked') + ' style=\"margin-right:6px;\">' +
                        node.label.split('\\n')[0] + '</label>';
        nodeListDiv.appendChild(div);
      });

      nodeListDiv.onchange = function(e) {
        if(e.target.classList.contains('node-check')) {
          var id = e.target.getAttribute('data-id'), n = nodesDS.get(id);
          n.hidden = !e.target.checked;
          nodesDS.update(n);
        }
      };

      visContent.querySelector('#sel_all').onclick = function() { 
        nodeListDiv.querySelectorAll('.node-check').forEach(c => { var n = nodesDS.get(c.getAttribute('data-id')); n.hidden = false; nodesDS.update(n); c.checked = true; });
      };
      visContent.querySelector('#sel_none').onclick = function() { 
        nodeListDiv.querySelectorAll('.node-check').forEach(c => { var n = nodesDS.get(c.getAttribute('data-id')); n.hidden = true; nodesDS.update(n); c.checked = false; });
      };

      var forceFit = function() { 
        var h = el.clientHeight;
        if(h > 0) {
          network.setSize('100%', h + 'px'); 
          network.redraw(); 
          network.fit(); 
        }
      };
      window.addEventListener('resize', forceFit);
      setTimeout(forceFit, 100); setTimeout(forceFit, 1500);
    }
    "

    g$x$interactive_options <- interactive_options
    g$x$min_weight    <- min_weight
    g$x$max_weight    <- max_w
    g$x$layouts       <- .get_all_layouts(wmatrix, vertex, area_attr)
    g$x$sim_data      <- sim_data
    g$x$initial_palette <- color
    g$x$initial_layout <- layout
    g$x$hide_direct    <- hide_direct
    g$x$color_palette_js <- list(
      "red/green" = c("#F52722", "#A5D610", "#999999", "#FFFF00"),
      "grey scale" = c("#808080", "#ffffff", "#f2f2f2", "#e5e5e5"),
      "colorblind" = c("#D55E00", "#0173B2", "#CC79A7", "#F0E442"),
      "pastel" = c("#f1677c", "#98FB98", "#F0F8FF", "#fcf087"),
      "dark" = c("#8B0000", "#006400", "#696969", "#DAA520"),
      "viridis" = c("#440154", "#35b779", "#31688e", "#fde725")
    )
    
    if (inherits(wimp, "wimp") && !is.null(wimp$vertices)) {
      w_verts <- wimp$vertices
      raw_cat_cols <- names(w_verts)[sapply(w_verts, function(x) is.character(x) || is.factor(x))]
      g$x$cat_cols <- setdiff(raw_cat_cols, c("left_pole", "right_pole", "label", "id", "self_pole", "ideal_pole"))
    } else {
      g$x$cat_cols <- c()
    }
    
    # Pass dictionary to JS
    g$x$dict <- t
    

    
    js_panel <- gsub("WIMP_EXPORT_NAME", export_name, js_panel)
    g <- g %>% htmlwidgets::onRender(js_panel)
  }

  return(g)
}


# Ideal Digraph -----------------------------------------------------------

#' Ideal digraph -- idealdigraph()
#'
#' @description Plot the ideal self based on the constructs and their relations.
#'
#' @param wimp Subject's WimpGrid. It must be a "wimp" S3 object
#'        imported by the \\code{\\link{importwimp}} function.
#' @param hide_direct If TRUE, hide direct relationships between nodes. Default is
#'            FALSE.
#' @param layout Layout for the digraph. Options: "circle", "rtcircle", "tree",
#'        "graphopt", "mds", "grid" or "areas". Default is "circle".
#' @param ... Additional arguments passed to \\code{\\link{digraph}}.
#'
#' @author Alejandro Sanfeliciano
#'
#' @return A visNetwork graph
#'
#' @export
#'
#' @examples
#' idealdigraph(example_wimp)

idealdigraph <- function(wimp, hide_direct = FALSE, layout = "circle", ...) {

  # Extract ideal vector from new-format vertices
  ideal_vector <- wimp$vertices$ideal
  plot <- digraph(wimp = wimp, hide_direct = hide_direct,
                  vertex_vector = ideal_vector, layout = layout, ...)
  return(plot)
}

# Simulation Digraph ---------------------------------------------------------

#' Simulation digraph -- simdigraph()
#'
#' @description Plot the hypothetical self for a given scenario iteration.
#'
#' @param scn A scenario matrix ("scn" S3 object) from
#'            \code{\link{scenariomatrix}}.
#' @param niter Iteration index to display (0-based). Default is 0.
#' @param ... Additional arguments passed to \code{\link{digraph}}.
#'
#' @author Alejandro Sanfeliciano
#'
#' @return A visNetwork graph
#'
#' @export
#'
#' @examples
#' scn <- scenariomatrix(example_wimp, rep(1,5))
#' simdigraph(scn, niter = 2)

simdigraph <- function(scn, niter = 0, ...) {
  # 1. Determine input type and build simulation metadata
  if (inherits(scn, "wimp")) {
    wimp <- .align_wimp(scn, exclude_dilemmatics = FALSE)
    sim_data <- list(
      initial_self    = wimp$vertices$self,
      initial_ideal   = wimp$vertices$ideal,
      act_vector      = rep(0, nrow(wimp$vertices)),
      weights         = as.matrix(wimp$global$weight_matrix),
      threshold       = "saturation",
      max_iter        = 10,
      lpoles          = wimp$vertices$left_pole,
      rpoles          = wimp$vertices$right_pole
    )
    vertex_vector <- wimp$vertices$self
    return(digraph(wimp = wimp, vertex_vector = vertex_vector, sim_data = sim_data, ...))
    
  } else if (inherits(scn, "scn")) {
    sim_data <- list(
      initial_self    = scn$self$self,
      initial_ideal   = scn$self$ideal,
      act_vector      = scn$params$act_vector,
      weights         = as.matrix(scn$weights),
      threshold       = scn$params$threshold,
      max_iter        = scn$params$max_iter,
      lpoles          = scn$constructs$left_pole,
      rpoles          = scn$constructs$right_pole
    )
    vertex_vector <- scn$values[niter + 1, ]
    return(digraph(wimp = scn, vertex_vector = vertex_vector, sim_data = sim_data, ...))
    
  } else {
    stop("Input must be a 'wimp' or 'scn' object.")
  }
}

# Weight Matrix Heatmap --------------------------------------------------------

#' Weight Matrix Heatmap -- weight_heatmap()
#'
#' @description Creates an interactive plotly heatmap of the weight matrix, 
#'              reoriented towards the ideal self poles.
#'
#' @param wimp Subject's WimpGrid object. It must be a "wimp" S3 object
#'        imported by the \code{\link{importwimp}} function.
#' @param palette Character. Color palette for the heatmap. Options: \code{"Redgreen"}
#'   (default, red-white-green), \code{"viridis"}, \code{"plasma"}, etc.
#'
#' @return A plotly heatmap object.
#'
#' @author Alejandro Sanfeliciano
#'
#' @import plotly
#' @export
#'
#' @examples
#' weight_heatmap(example_wimp)
#'
weight_heatmap <- function(wimp, palette = "Redgreen", lang = "en") {
  
  # 1. Extract data
  if (!inherits(wimp, "wimp")) {
    stop("'wimp' must be an object of class 'wimp'.")
  }
  
  wmatrix <- as.matrix(wimp$global$weight_matrix)
  ideal_vector <- wimp$vertices$ideal
  l_poles <- wimp$vertices$left_pole
  r_poles <- wimp$vertices$right_pole
  
  # 2. Reorient matrix
  # Use signs of ideal vector to reorient relations towards the "desired" direction
  signs <- sign(ideal_vector)
  signs[signs == 0] <- 1
  
  # Structural reorientation: diag(signs) %*% W %*% diag(signs)
  wmatrix_ideal <- diag(signs, nrow = length(signs)) %*% wmatrix %*% diag(signs, nrow = length(signs))
  
  # 3. Prepare labels (pole towards ideal)
  labels <- sapply(seq_along(ideal_vector), function(i) {
    if (ideal_vector[i] > 0) {
      r_poles[i]
    } else if (ideal_vector[i] < 0) {
      l_poles[i]
    } else {
      r_poles[i]
    }
  })
  
  # 4. Palette management
  if (palette == "Redgreen") {
    palette <- list(c(0, "#F52722"), c(0.5, "white"), c(1, "#A5D610"))
  }
  
  # 5. Create Heatmap
  p <- plot_ly(
    x = labels,
    y = labels,
    z = wmatrix_ideal,
    type = "heatmap",
    colorscale = palette,
    zmin = -2,
    zmax = 2,
    xgap = 0,
    ygap = 0,
    hovertemplate = if (lang == "es") {
      paste(
        "<b>Origen:</b> %{y}<br>",
        "<b>Destino:</b> %{x}<br>",
        "<b>Peso:</b> %{z:.3f}<extra></extra>"
      )
    } else {
      paste(
        "<b>From:</b> %{y}<br>",
        "<b>To:</b> %{x}<br>",
        "<b>Weight:</b> %{z:.3f}<extra></extra>"
      )
    },
    colorbar = list(title = "", len = 0.75, y = 0.45, yanchor = "middle")
  )
  
  # 6. Add dotted grid lines and dilemmatic highlighting
  n_con <- length(labels)
  shapes <- list()
  is_dilemmatic <- ideal_vector == 0
  
  # Color/style for dilemmatic highlighting
  d_line_color <- "#FFC107" # Orange-Yellow
  d_fill_color <- "rgba(255, 193, 7, 0.10)" # More subtle transparent tint
  
  # A. Add Tinted Rectangles for Dilemmatic rows/columns
  for (i in seq_along(is_dilemmatic)) {
    if (is_dilemmatic[i]) {
      # Column tint
      shapes[[length(shapes) + 1]] <- list(
        type = "rect", x0 = i - 1.5, x1 = i - 0.5, y0 = -0.5, y1 = n_con - 0.5,
        fillcolor = d_fill_color, line = list(width = 0), layer = "above"
      )
      # Row tint
      shapes[[length(shapes) + 1]] <- list(
        type = "rect", x0 = -0.5, x1 = n_con - 0.5, y0 = i - 1.5, y1 = i - 0.5,
        fillcolor = d_fill_color, line = list(width = 0), layer = "above"
      )
    }
  }

  # B. Add dotted grid lines (between cells)
  if (n_con > 1) {
    for (i in 1:(n_con - 1)) {
      # Highlight vertical line if either adjacent construct is dilemmatic
      use_strong <- is_dilemmatic[i] | is_dilemmatic[i+1]
      l_color <- if (use_strong) d_line_color else "rgba(0,0,0,0.15)"
      l_width <- if (use_strong) 2 else 1
      
      shapes[[length(shapes) + 1]] <- list(
        type = "line", x0 = i - 0.5, x1 = i - 0.5, y0 = -0.5, y1 = n_con - 0.5,
        line = list(color = l_color, width = l_width, dash = "dot"), layer = "above"
      )
    }
    for (i in 1:(n_con - 1)) {
      # Highlight horizontal line if either adjacent construct is dilemmatic
      use_strong <- is_dilemmatic[i] | is_dilemmatic[i+1]
      l_color <- if (use_strong) d_line_color else "rgba(0,0,0,0.15)"
      l_width <- if (use_strong) 2 else 1
      
      shapes[[length(shapes) + 1]] <- list(
        type = "line", x0 = -0.5, x1 = n_con - 0.5, y0 = i - 0.5, y1 = i - 0.5,
        line = list(color = l_color, width = l_width, dash = "dot"), layer = "above"
      )
    }
  }

  p <- p %>%
  layout(
    title = "",
    xaxis = list(
      title = list(text = if (lang == "es") "<b>Efecto sobre (Destino)</b>" else "<b>Effect on (To)</b>", font = list(size = 14), standoff = 25),
      tickangle = -45,
      tickfont = list(size = 10),
      showline = TRUE, mirror = TRUE, linecolor = "black", linewidth = 1,
      showgrid = FALSE, zeroline = FALSE
    ),
    yaxis = list(
      title = list(text = if (lang == "es") "<b>Influencia de (Origen)</b>" else "<b>Influence of (From)</b>", font = list(size = 14), standoff = 25),
      autorange = "reversed",
      tickfont = list(size = 10),
      showline = TRUE, mirror = TRUE, linecolor = "black", linewidth = 1,
      showgrid = FALSE, zeroline = FALSE
    ),
    shapes = shapes,
    annotations = list(
      list(
        x = 1.05, y = 1.0,
        text = paste0("<b>\u03c1(G) = ", round(density_index(wimp), 3), "</b>"),
        showarrow = FALSE,
        xref = "paper", yref = "paper",
        xanchor = "left", yanchor = "top",
        font = list(size = 12)
      )
    ),
    margin = list(l = 120, r = 120, b = 120, t = 60)
  ) %>%
  config(displayModeBar = FALSE)
  
  hm_data <- list(
    dict = wt_i18n(lang),
    orig_labels = labels,
    orig_matrix = wmatrix_ideal,
    initial_palette = palette,
    is_dilemmatic = is_dilemmatic
  )
  
  js_hm_panel <- "
    function(el, p_x, data) {
      var x = data;
      // Create Settings Modal
      var settingsModal = document.createElement('div');
      settingsModal.id = 'hm_settings_modal';
      Object.assign(settingsModal.style, {
        position: 'absolute', top: '10px', right: '10px',
        width: '90%', maxWidth: '300px', maxHeight: '80vh', overflowY: 'auto', boxSizing: 'border-box',
        backgroundColor: '#fff', zIndex: '2000', padding: '20px', borderRadius: '8px', 
        boxShadow: '0 4px 20px rgba(0,0,0,0.2)', border: '1px solid #eaeaea', display: 'none', 
        fontFamily: 'Inter, Roboto, sans-serif'
      });
      
      var hmHTML = '<div style=\"display:flex; justify-content:space-between; align-items:center; border-bottom:1px solid #eaeaea; padding-bottom:10px; margin-bottom:15px;\">' +
                   '<h3 style=\"margin:0; color:#444; font-size:14px;\">' + x.dict.hm_settings + '</h3>' +
                   '<span id=\"close_hm_settings\" style=\"cursor:pointer; font-size:20px; font-weight:bold; color:#888; line-height:1;\">&times;</span>' +
                   '</div>';
                   
      hmHTML += '<div style=\"margin-bottom:15px;\">' +
                '<label style=\"display:block; margin-bottom:5px; font-weight:bold; color:#444; font-size:13px;\">' + x.dict.color_palette + '</label>' +
                '<select id=\"hm_palette_sel\" style=\"width:100%; padding:4px; border-radius:4px;\">' +
                '<option value=\"Redgreen\"' + (x.initial_palette === \"Redgreen\" ? ' selected' : '') + '>' + x.dict.pal_redgreen + '</option>' +
                '<option value=\"Redblue\">' + x.dict.pal_redblue + '</option>' +
                '<option value=\"Orangepurple\">' + x.dict.pal_orangepurple + '</option>' +
                '<option value=\"Greyscale\">' + x.dict.pal_greyscale + '</option>' +
                '</select></div>';
                
      hmHTML += '<div style=\"margin-bottom:15px;\">' +
                '<label style=\"display:block; margin-bottom:5px; font-weight:bold; color:#444; font-size:13px;\">' + x.dict.sort_by + '</label>' +
                '<select id=\"hm_sort_sel\" style=\"width:100%; padding:4px; border-radius:4px;\">' +
                '<option value=\"original\" selected>' + x.dict.sort_original + '</option>' +
                '<option value=\"weight\">' + x.dict.sort_weight + '</option>' +
                '<option value=\"connectivity\">' + x.dict.sort_connect + '</option>' +
                '</select></div>';
                
      hmHTML += '<div style=\"margin-bottom:15px;\">' +
                '<label style=\"display:block; margin-bottom:5px; font-weight:bold; color:#444; font-size:13px;\">' + x.dict.filter_constructs + '</label>' +
                '<div id=\"hm_filter_list\" style=\"max-height:120px; overflow-y:auto; border:1px solid #ddd; padding:5px; border-radius:4px; font-size:12px; background:#f9f9f9;\"></div>' +
                '</div>';
                
      hmHTML += '<div style=\"margin-bottom:15px; border-top:1px solid #f0f0f0; padding-top:10px;\">' +
                '<label style=\"display:flex; align-items:center; font-size:12px; cursor:pointer; color:#555; margin-bottom:8px;\">' +
                '<input type=\"checkbox\" id=\"hm_density_check\" checked style=\"margin-right:6px;\"> ' + x.dict.show_density + '</label>' +
                '<label style=\"display:flex; align-items:center; font-size:12px; cursor:pointer; color:#555;\">' +
                '<input type=\"checkbox\" id=\"hm_values_check\" style=\"margin-right:6px;\"> ' + x.dict.show_values + '</label>' +
                '</div>';
                
      settingsModal.innerHTML = hmHTML;
      var container = el.closest('.wt-tab-content') || el; container.appendChild(settingsModal);
      
      var filterContainer = settingsModal.querySelector('#hm_filter_list');
      x.orig_labels.forEach(function(lbl, idx) {
         var div = document.createElement('div');
         div.style.marginBottom = '4px';
         div.innerHTML = '<label style=\"cursor:pointer; display:flex; align-items:center; color:#555;\"><input type=\"checkbox\" checked value=\"' + idx + '\" class=\"hm-construct-cb\" style=\"margin-right:6px;\"> ' + lbl + '</label>';
         filterContainer.appendChild(div);
      });
      
      settingsModal.querySelector('#close_hm_settings').onclick = function() { settingsModal.style.display = 'none'; };
      
      // Settings Button
      var settingsBtn = document.createElement('div');
      
      var hmBtnContainer = el.closest('.wt-tab-content') ? el.closest('.wt-tab-content').querySelector('#hm_btn_container') : null;
      
      if (hmBtnContainer) {
        Object.assign(settingsBtn.style, {
          backgroundColor: 'rgba(255, 255, 255, 0.95)', width: '32px', height: '32px',
          borderRadius: '6px', boxShadow: '0 2px 10px rgba(0,0,0,0.1)', border: '1px solid #ddd',
          display: 'flex', alignItems: 'center', justifyContent: 'center', cursor: 'pointer', transition: 'all 0.2s'
        });
      } else {
        Object.assign(settingsBtn.style, {
          position: 'absolute', bottom: '15px', right: '15px', zIndex: '1000',
          backgroundColor: 'rgba(255, 255, 255, 0.95)', width: '32px', height: '32px',
          borderRadius: '6px', boxShadow: '0 2px 10px rgba(0,0,0,0.1)', border: '1px solid #ddd',
          display: 'flex', alignItems: 'center', justifyContent: 'center', cursor: 'pointer', transition: 'all 0.2s'
        });
      }
      
      settingsBtn.innerHTML = \"<svg width='60%' height='60%' viewBox='0 0 24 24' fill='none' stroke='#333' stroke-width='2.5' stroke-linecap='round' stroke-linejoin='round'><circle cx='12' cy='12' r='3'></circle><path d='M19.4 15a1.65 1.65 0 0 0 .33 1.82l.06.06a2 2 0 0 1 0 2.83 2 2 0 0 1-2.83 0l-.06-.06a1.65 1.65 0 0 0-1.82-.33 1.65 1.65 0 0 0-1 1.51V21a2 2 0 0 1-2 2 2 2 0 0 1-2-2v-.09A1.65 1.65 0 0 0 9 19.4a1.65 1.65 0 0 0-1.82.33l-.06.06a2 2 0 0 1-2.83 0 2 2 0 0 1 0-2.83l.06-.06a1.65 1.65 0 0 0 .33-1.82 1.65 1.65 0 0 0-1.51-1H3a2 2 0 0 1-2-2 2 2 0 0 1 2-2h.09A1.65 1.65 0 0 0 4.6 9a1.65 1.65 0 0 0-.33-1.82l-.06-.06a2 2 0 0 1 0-2.83 2 2 0 0 1 2.83 0l.06.06a1.65 1.65 0 0 0 1.82.33H9a1.65 1.65 0 0 0 1-1.51V3a2 2 0 0 1 2-2 2 2 0 0 1 2 2v.09a1.65 1.65 0 0 0 1 1.51 1.65 1.65 0 0 0 1.82-.33l.06-.06a2 2 0 0 1 2.83 0 2 2 0 0 1 0 2.83l-.06.06a1.65 1.65 0 0 0-.33 1.82V9a1.65 1.65 0 0 0 1.51 1H21a2 2 0 0 1 2 2 2 2 0 0 1-2 2h-.09a1.65 1.65 0 0 0-1.51 1z'></path></svg>\";
      settingsBtn.title = \"Ajustes\";
      settingsBtn.onmouseover = function() { this.style.backgroundColor = '#f5f5f5'; };
      settingsBtn.onmouseout = function() { this.style.backgroundColor = 'rgba(255, 255, 255, 0.95)'; };
      settingsBtn.onclick = function() { settingsModal.style.display = (settingsModal.style.display === 'block' ? 'none' : 'block'); };
      
      if (hmBtnContainer) {
        var placeholder = hmBtnContainer.querySelector('#hm_settings_placeholder');
        if (placeholder) {
          hmBtnContainer.replaceChild(settingsBtn, placeholder);
        } else {
          hmBtnContainer.appendChild(settingsBtn);
        }
      } else {
        var container = el.closest('.wt-tab-content') || el;
        container.appendChild(settingsBtn);
      }

      // Save original annotations to toggle density
      var origAnnotations = el.layout.annotations ? JSON.parse(JSON.stringify(el.layout.annotations)) : [];

      // Logic functions
      var updateHeatmap = function() {
        var pal = settingsModal.querySelector('#hm_palette_sel').value;
        var srt = settingsModal.querySelector('#hm_sort_sel').value;
        var showDens = settingsModal.querySelector('#hm_density_check').checked;
        var showVals = settingsModal.querySelector('#hm_values_check').checked;
        
        var cs = [[0, '#F52722'], [0.5, 'white'], [1, '#A5D610']]; // Default Redgreen
        if(pal === 'Redblue') cs = [[0, '#F52722'], [0.5, 'white'], [1, '#2272F5']];
        if(pal === 'Orangepurple') cs = [[0, '#F58222'], [0.5, 'white'], [1, '#9B22F5']];
        if(pal === 'Greyscale') cs = [[0, '#000000'], [0.5, 'white'], [1, '#000000']]; // V-shaped gradient for absolute values
        
        var activeIndices = [];
        settingsModal.querySelectorAll('.hm-construct-cb').forEach(function(cb) {
           if(cb.checked) activeIndices.push(parseInt(cb.value));
        });
        
        // Sorting logic
        var labels = x.orig_labels.slice();
        var full_n = labels.length;
        var indices = activeIndices.slice();
        var n = indices.length;
        
        if (srt === 'weight') {
           var scores = indices.map(i => {
              var sum = 0;
              for(var j=0; j<full_n; j++) sum += Math.abs(x.orig_matrix[i][j]) + Math.abs(x.orig_matrix[j][i]);
              return sum;
           });
           indices.sort((a,b) => scores[b] - scores[a]);
        } else if (srt === 'connectivity') {
           var scores = indices.map(i => {
              var count = 0;
              for(var j=0; j<full_n; j++) {
                 if(Math.abs(x.orig_matrix[i][j]) > 0.001) count++;
                 if(Math.abs(x.orig_matrix[j][i]) > 0.001) count++;
              }
              return count;
           });
           indices.sort((a,b) => scores[b] - scores[a]);
        }
        
        var new_labels = indices.map(i => labels[i]);
        var new_z = [];
        var new_text = [];
        for(var i=0; i<n; i++) {
           var row_z = [];
           var row_txt = [];
           for(var j=0; j<n; j++) {
               var val = x.orig_matrix[indices[i]][indices[j]];
               row_z.push(val);
               if (pal === 'Greyscale' && val < -0.01) {
                   row_txt.push('-'); // Minus indicator
               } else {
                   row_txt.push('');
               }
           }
           new_z.push(row_z);
           new_text.push(row_txt);
        }
        
        var restyleData = {
           colorscale: [cs],
           x: [new_labels],
           y: [new_labels],
           z: [new_z],
           text: [new_text]
        };
        
        if (showVals) {
           restyleData.texttemplate = ['%{z:.2f}'];
           restyleData.textfont = [{color: 'black', size: 10}];
        } else if (pal === 'Greyscale') {
           restyleData.texttemplate = ['%{text}'];
           restyleData.textfont = [{color: '#F52722', size: 18, family: 'Arial'}];
        } else {
           restyleData.texttemplate = [null];
        }
        
        Plotly.restyle(el, restyleData, [0]);
        var new_shapes = [];
        var d_line_color = '#FFC107';
        var d_fill_color = 'rgba(255, 193, 7, 0.10)';
        
        for (var i = 0; i < n; i++) {
            if (x.is_dilemmatic[indices[i]]) {
                new_shapes.push({
                   type: 'rect', x0: i - 0.5, x1: i + 0.5, y0: -0.5, y1: n - 0.5,
                   fillcolor: d_fill_color, line: {width: 0}, layer: 'above'
                });
                new_shapes.push({
                   type: 'rect', x0: -0.5, x1: n - 0.5, y0: i - 0.5, y1: i + 0.5,
                   fillcolor: d_fill_color, line: {width: 0}, layer: 'above'
                });
            }
        }
        
        if (n > 1) {
            for (var i = 0; i < n - 1; i++) {
                var use_strong = x.is_dilemmatic[indices[i]] || x.is_dilemmatic[indices[i+1]];
                var l_color = use_strong ? d_line_color : 'rgba(0,0,0,0.15)';
                var l_width = use_strong ? 2 : 1;
                
                new_shapes.push({
                   type: 'line', x0: i + 0.5, x1: i + 0.5, y0: -0.5, y1: n - 0.5,
                   line: {color: l_color, width: l_width, dash: 'dot'}, layer: 'above'
                });
                new_shapes.push({
                   type: 'line', x0: -0.5, x1: n - 0.5, y0: i + 0.5, y1: i + 0.5,
                   line: {color: l_color, width: l_width, dash: 'dot'}, layer: 'above'
                });
            }
        }

        Plotly.relayout(el, {
           annotations: showDens ? origAnnotations : [],
           shapes: new_shapes,
           'xaxis.showticklabels': true,
           'yaxis.showticklabels': true
        });
      };
      
      settingsModal.querySelector('#hm_palette_sel').onchange = updateHeatmap;
      settingsModal.querySelector('#hm_sort_sel').onchange = updateHeatmap;
      settingsModal.querySelector('#hm_density_check').onchange = updateHeatmap;
      settingsModal.querySelector('#hm_values_check').onchange = updateHeatmap;
      settingsModal.querySelectorAll('.hm-construct-cb').forEach(function(cb) {
         cb.onchange = updateHeatmap;
      });
    }
"
  
  p <- p %>% htmlwidgets::onRender(js_hm_panel, data = hm_data)
  
  return(p)
}
