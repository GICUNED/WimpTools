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
                    interactive_options = TRUE, sim_data = NULL, ...) {

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
      return(list(color = character(0), dashes = logical(0)))
    }

    if (color != "grey scale") {
      list(
        color = sapply(weight_vector, function(x) {
          ifelse(x > 0, "grey", "#CD5C5C")
        }),
        dashes = rep(FALSE, length(weight_vector))
      )
    } else {
      list(
        color = rep("grey", length(weight_vector)),
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
  if (hide_direct) {
    logical_dilemmatic <- ideal_vector == 0
    wmatrix[wmatrix > 0] <- 0
    wmatrix[logical_dilemmatic, ] <- 0
    wmatrix[, logical_dilemmatic] <- 0
  }

  # Extract all potential edges and their bidirectional status
  all_edges_raw <- .extract_edges(wmatrix)
  all_edge_curved <- .detect_bidirectional_edges(wmatrix)
  
  max_w <- (if (nrow(all_edges_raw) > 0) max(abs(all_edges_raw$weight)) else 1) + 0.01

  edges_raw <- all_edges_raw
  edge_curved <- all_edge_curved

  if (nrow(edges_raw) == 0) {
    edges <- data.frame(
      from = integer(0), to = integer(0), width = numeric(0),
      arrows = character(0), dashes = logical(0), smooth = logical(0),
      color = character(0), title = numeric(0), weight = numeric(0),
      hidden = logical(0),
      stringsAsFactors = FALSE
    )
  } else {
    edge_props <- .calculate_edge_properties(edges_raw$weight, color)
    
    # Use subtle curves for bidirectional edges
    smooth_list <- lapply(edge_curved, function(sc) {
      if (sc) list(enabled = TRUE, type = "curvedCW", roundness = 0.15) else list(enabled = FALSE)
    })

    edges <- data.frame(
      from = as.character(edges_raw$from),
      to = as.character(edges_raw$to),
      width = 2 * abs(edges_raw$weight),
      arrows = "to",
      dashes = edge_props$dashes,
      smooth = smooth_list,
      color = edge_props$color,
      title = round(edges_raw$weight, 2),
      weight = edges_raw$weight,
      orig_dashes = edge_props$dashes,
      orig_color = edge_props$color,
      hidden = if (interactive_options) abs(edges_raw$weight) < min_weight else FALSE,
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
      visOptions(manipulation = list(enabled = TRUE, addNode = FALSE, addEdge = FALSE, 
                                     editNode = FALSE, editEdge = FALSE, 
                                     deleteNode = TRUE, deleteEdge = TRUE),
                 highlightNearest = list(enabled = TRUE, degree = 0,
                                         labelOnly = TRUE),
                 selectedBy = list(variable = "group", main = "All")) %>%
      visInteraction(navigationButtons = FALSE, multiselect = TRUE) %>%
      visPhysics(enabled = FALSE) %>%
      visNodes(font = list(align = "center", multi = TRUE, vadjust = 0))
  } else if (layout == "rtcircle") {
    # Circular tree layout
    g <- visNetwork(vertex, edges, height = height, width = width) %>%
      visIgraphLayout(layout = "layout_as_tree", circular = TRUE) %>%
      visOptions(manipulation = list(enabled = TRUE, addNode = FALSE, addEdge = FALSE, 
                                     editNode = FALSE, editEdge = FALSE, 
                                     deleteNode = TRUE, deleteEdge = TRUE),
                 highlightNearest = list(enabled = TRUE, degree = 0,
                                         labelOnly = TRUE),
                 selectedBy = list(variable = "group", main = "All")) %>%
      visInteraction(navigationButtons = FALSE, multiselect = TRUE) %>%
      visNodes(font = list(align = "center", multi = TRUE, vadjust = 0))
  } else {
    # Standard igraph layouts
    g <- visNetwork(vertex, edges, height = height, width = width) %>%
      visIgraphLayout(layout = layout_name, randomSeed = 33) %>%
      visOptions(manipulation = list(enabled = TRUE, addNode = FALSE, addEdge = FALSE, 
                                     editNode = FALSE, editEdge = FALSE, 
                                     deleteNode = TRUE, deleteEdge = TRUE),
                 highlightNearest = list(enabled = TRUE, degree = 0,
                                         labelOnly = TRUE),
                 selectedBy = list(variable = "group", main = "All")) %>%
      visInteraction(navigationButtons = FALSE, multiselect = TRUE) %>%
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
      el.style.height = '90vh';
      el.style.minHeight = '600px';
      var network = this.network;
      var container = el;
      
      var currentDistMult = 1.0;
      var currentSizeMult = 1.0;
      var currentTextSize = 20;

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
        var currentIdx = (window._simCurrentI !== undefined) ? window._simCurrentI : 0;
        var vals = simVals || (window._simHistory && window._simHistory[currentIdx]);
        var nodesDS = network.body.data.nodes;
        
        var updates = nodesDS.get().map(function(node, i) {
          var val = (vals && vals[i] !== undefined) ? vals[i] : (node.self_val || 0);
          var ideal = (x.sim_data) ? x.sim_data.initial_ideal[i] : node.ideal_val;
          var c = getPaletteColor(val, ideal, scheme);
          var label = node.label;
          if (x.sim_data && window._simHistory) {
            label = (val < 0) ? x.sim_data.lpoles[i] : (val > 0 ? x.sim_data.rpoles[i] : x.sim_data.lpoles[i] + ' - ' + x.sim_data.rpoles[i]);
          }
          var baseSize = (x.sim_data && window._simHistory) ? (20 + (30 * Math.abs(val))) : (node.raw_size || 20);
          var finalSize = baseSize * currentSizeMult;
          var vadjust = - (finalSize * 1.1 + currentTextSize * 0.8);
          return {
            id: node.id,
            color: { background: c, border: darkenColor(c, 0.4), highlight: { background: c, border: darkenColor(c, 0.6) } },
            size: finalSize, label: label, shape: 'dot',
            font: { vadjust: vadjust, size: currentTextSize, face: 'Segoe UI', color: '#000000', strokeWidth: 3, strokeColor: '#ffffff' }
          };
        });
        nodesDS.update(updates);
        
        var edgesDS = network.body.data.edges;
        var scheme2 = visContent.querySelector('#palette_sel').value;
        var edgeUpdates = edgesDS.get().map(function(edge) {
          if(scheme2 === 'grey scale') return {id: edge.id, color: '#999999', dashes: edge.weight < 0};
          return {id: edge.id, color: edge.orig_color, dashes: edge.orig_dashes};
        });
        edgesDS.update(edgeUpdates);
      };

      // ── Edge flow flash: color edges by activation contribution ──────────
      var flashEdgeFlow = function(fromIdx, toIdx) {
        if(!window._simHistory || !window._simHistory[fromIdx] || !window._simHistory[toIdx]) return;
        var prevVals = window._simHistory[fromIdx];
        var nextVals = window._simHistory[toIdx];
        var deltas   = nextVals.map(function(v, i) { return v - prevVals[i]; });
        var weights  = x.sim_data.weights; // weights[from][to]
        var THRESHOLD = 0.04;
        var edgesDS  = network.body.data.edges;
        var flowUpdates = edgesDS.get().map(function(edge) {
          // edge.from / edge.to are node IDs (strings); find their index
          var nodesDS = network.body.data.nodes;
          var allNodes = nodesDS.get();
          var srcIdx = allNodes.findIndex(function(n) { return n.id === edge.from; });
          var dstIdx = allNodes.findIndex(function(n) { return n.id === edge.to;   });
          if(srcIdx < 0 || dstIdx < 0) return {id: edge.id};
          var flow = (weights[srcIdx] && weights[srcIdx][dstIdx] !== undefined)
            ? weights[srcIdx][dstIdx] * deltas[srcIdx]
            : 0;
          var color;
          if     (flow >  THRESHOLD) color = {color: '#4CAF50', highlight: '#4CAF50', hover: '#4CAF50'};
          else if(flow < -THRESHOLD) color = {color: '#E53935', highlight: '#E53935', hover: '#E53935'};
          else                       color = edge.orig_color;
          return {id: edge.id, color: color};
        });
        edgesDS.update(flowUpdates);
      };

      // ── Smooth tween between two iteration states ─────────────────────────
      var _tweenRAF = null;
      var tweenToIteration = function(fromIdx, toIdx, durationMs) {
        if(_tweenRAF) { cancelAnimationFrame(_tweenRAF); _tweenRAF = null; }
        if(!window._simHistory || !window._simHistory[fromIdx] || !window._simHistory[toIdx]) {
          window._simCurrentI = toIdx;
          refreshNodes();
          return;
        }
        var prevVals = window._simHistory[fromIdx].slice();
        var nextVals = window._simHistory[toIdx];
        var nodesDS  = network.body.data.nodes;
        var allNodes = nodesDS.get();
        var scheme   = visContent.querySelector('#palette_sel').value;
        var t0 = null;
        // Flash edge flow at start of tween
        flashEdgeFlow(fromIdx, toIdx);
        var step = function(ts) {
          if(!t0) t0 = ts;
          var t = Math.min((ts - t0) / durationMs, 1);
          // Ease-in-out cubic
          var ease = t < 0.5 ? 4*t*t*t : 1 - Math.pow(-2*t+2, 3)/2;
          // Interpolate node values and render
          var interpVals = nextVals.map(function(nv, i) { return prevVals[i] + ease * (nv - prevVals[i]); });
          var updates = allNodes.map(function(node, i) {
            var val   = interpVals[i] !== undefined ? interpVals[i] : (node.self_val || 0);
            var ideal = x.sim_data ? x.sim_data.initial_ideal[i] : node.ideal_val;
            var c     = getPaletteColor(val, ideal, scheme);
            var baseSize = 20 + 30 * Math.abs(val);
            var finalSize = baseSize * currentSizeMult;
            var vadjust   = -(finalSize * 1.1 + currentTextSize * 0.8);
            var label = (val < 0) ? x.sim_data.lpoles[i] : (val > 0 ? x.sim_data.rpoles[i] : x.sim_data.lpoles[i] + ' - ' + x.sim_data.rpoles[i]);
            return {
              id: node.id,
              color: { background: c, border: darkenColor(c, 0.4), highlight: { background: c, border: darkenColor(c, 0.6) } },
              size: finalSize, label: label, shape: 'dot',
              font: { vadjust: vadjust, size: currentTextSize, face: 'Segoe UI', color: '#000000', strokeWidth: 3, strokeColor: '#ffffff' }
            };
          });
          nodesDS.update(updates);
          if(t < 1) {
            _tweenRAF = requestAnimationFrame(step);
          } else {
            _tweenRAF = null;
            window._simCurrentI = toIdx;
            refreshNodes(); // snap to exact final state + restore edge colors
          }
        };
        _tweenRAF = requestAnimationFrame(step);
      };


      // --- Export Logic ---
      var exportPNG = function() {
        var canvas = container.getElementsByTagName('canvas')[0];
        if(!canvas) return;
        var link = document.createElement('a');
        link.download = 'WimpTools_Digraph_' + new Date().getTime() + '.png';
        link.href = canvas.toDataURL('image/png', 1.0);
        link.click();
      };

      // --- Helper to create panels ---
      var createPanel = function(id, title, positionStyles) {
        var p = document.createElement('div');
        p.id = id;
        Object.assign(p.style, {
          position: 'absolute', zIndex: '1000', backgroundColor: 'rgba(255, 255, 255, 0.95)',
          padding: '10px', borderRadius: '8px', boxShadow: '0 2px 15px rgba(0,0,0,0.15)',
          border: '1px solid #ddd', fontFamily: 'Segoe UI, Tahoma, sans-serif', fontSize: '12px',
          width: '220px', maxHeight: '40px', overflowY: 'hidden', transition: 'all 0.3s ease'
        }, positionStyles);

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
        toggleIcon.className = 'toggle-icon';
        toggleIcon.innerText = '+';
        header.appendChild(toggleIcon);
        p.appendChild(header);

        var content = document.createElement('div');
        content.style.display = 'none';
        content.style.marginTop = '10px';
        content.style.borderTop = '1px solid #eee';
        content.style.paddingTop = '10px';
        p.appendChild(content);

        header.onclick = function() {
          var isHidden = content.style.display === 'none';
          content.style.display = isHidden ? 'block' : 'none';
          toggleIcon.innerText = isHidden ? '−' : '+';
          p.style.maxHeight = isHidden ? '90%' : '40px';
          p.style.overflowY = isHidden ? 'auto' : 'hidden';
        };

        container.appendChild(p);
        return content;
      };

      // --- Panels Initialization ---
      var visContent = createPanel('vis_panel', 'Visualization Options', {top: '10px', right: '10px'});
      
      // Minimalist Export Panel
      var exportPanel = document.createElement('div');
      Object.assign(exportPanel.style, {
        position: 'absolute', top: '10px', right: '265px', zIndex: '1000',
        backgroundColor: 'rgba(255, 255, 255, 0.95)', width: '32px', height: '32px',
        borderRadius: '6px', boxShadow: '0 2px 10px rgba(0,0,0,0.1)',
        border: '1px solid #ddd', display: 'flex', alignItems: 'center',
        justifyContent: 'center', cursor: 'pointer', transition: 'all 0.2s'
      });
      exportPanel.innerHTML = \"<svg width='18' height='18' viewBox='0 0 24 24' fill='none' stroke='#333' stroke-width='2.5' stroke-linecap='round' stroke-linejoin='round'><path d='M21 15v4a2 2 0 0 1-2 2H5a2 2 0 0 1-2-2v-4'></path><polyline points='7 10 12 15 17 10'></polyline><line x1='12' y1='15' x2='12' y2='3'></line></svg>\";
      exportPanel.title = 'Export PNG';
      exportPanel.onmouseover = function() { this.style.backgroundColor = '#f5f5f5'; };
      exportPanel.onmouseout = function() { this.style.backgroundColor = 'rgba(255, 255, 255, 0.95)'; };
      exportPanel.onclick = exportPNG;
      container.appendChild(exportPanel);

      // --- Visualization Content ---
      var visHTML = '<div style=\"margin-bottom:15px;\">' +
                    '<label style=\"display:block; margin-bottom:5px; font-weight:bold; color:#444;\">Color Palette</label>' +
                    '<select id=\"palette_sel\" style=\"width:100%; padding:4px; border-radius:4px;\">' +
                    Object.keys(x.color_palette_js).map(k => '<option value=\"' + k + '\"' + (k === x.initial_palette ? ' selected' : '') + '>' + k.charAt(0).toUpperCase() + k.slice(1) + '</option>').join('') +
                    '</select></div>';
      
      visHTML += '<div style=\"margin-bottom:15px;\">' +
                 '<label style=\"display:block; margin-bottom:5px; font-weight:bold; color:#444;\">Layout Algorithm</label>' +
                 '<select id=\"layout_sel\" style=\"width:100%; padding:4px; border-radius:4px;\">' +
                 Object.keys(x.layouts).map(k => '<option value=\"' + k + '\"' + (k === x.initial_layout ? ' selected' : '') + '>' + k.charAt(0).toUpperCase() + k.slice(1) + '</option>').join('') +
                 '</select></div>';
      
      visHTML += '<div style=\"margin-bottom:15px; border-top:1px solid #eee; padding-top:10px;\">' +
                 '<div style=\"display:flex; justify-content:space-between; margin-bottom:5px;\"><b style=\"color:#444;\">Edge Filter</b><span id=\"weight_val_txt\" style=\"font-family:monospace;\">0.00</span></div>' +
                 '<input type=\"range\" id=\"weight_slider\" min=\"0\" max=\"' + x.max_weight + '\" step=\"0.01\" value=\"0\" style=\"width:100%;\">' +
                 '</div>';

      visHTML += '<div style=\"margin-bottom:15px; border-top:1px solid #f0f0f0; padding-top:10px;\">' +
                 '<div style=\"display:flex; justify-content:space-between; margin-bottom:5px;\"><b style=\"color:#444;\">Node Spacing</b></div>' +
                 '<input type=\"range\" id=\"dist_slider\" min=\"0.5\" max=\"3\" step=\"0.1\" value=\"1\" style=\"width:100%;\">' +
                 '</div>';

      visHTML += '<div style=\"margin-bottom:15px;\">' +
                 '<div style=\"display:flex; justify-content:space-between; margin-bottom:5px;\"><b style=\"color:#444;\">Node Size</b></div>' +
                 '<input type=\"range\" id=\"size_slider\" min=\"0.5\" max=\"3\" step=\"0.1\" value=\"1\" style=\"width:100%;\">' +
                 '</div>';

      visHTML += '<div style=\"margin-bottom:15px;\">' +
                 '<div style=\"display:flex; justify-content:space-between; margin-bottom:5px;\"><b style=\"color:#444;\">Text Size</b></div>' +
                 '<input type=\"range\" id=\"text_size_slider\" min=\"10\" max=\"40\" step=\"1\" value=\"20\" style=\"width:100%;\">' +
                 '</div>';

      visHTML += '<div style=\"border-top:1px solid #eee; padding-top:10px; margin-bottom:15px;\">' +
                 '<button id=\"btn_reset\" style=\"width:100%; padding:6px; background:#f8f9fa; border:1px solid #ccc; border-radius:4px; font-weight:bold; color:#555; cursor:pointer;\">Reset Configuration</button>' +
                 '</div>';

      visHTML += '<div style=\"border-top:1px solid #eee; padding-top:10px;\">' +
                 '<b style=\"color:#444; display:block; margin-bottom:8px;\">Visible Constructs</b>' +
                 '<div style=\"display:flex; gap:5px; margin-bottom:8px;\">' +
                 '<button id=\"sel_all\" style=\"flex:1; font-size:10px; cursor:pointer;\">All</button>' +
                 '<button id=\"sel_none\" style=\"flex:1; font-size:10px; cursor:pointer;\">None</button></div>' +
                 '<div id=\"node_list\" style=\"max-height:150px; overflow-y:auto; border:1px solid #f0f0f0; padding:5px;\"></div></div>';
      
      visContent.innerHTML = visHTML;

      // --- Simulation Panel ---
      if (x.sim_data) {
        var sim = x.sim_data;
        window._simHistory = [];
        window._simCurrentI = 0;
        var targetSelf = [...sim.initial_self];
        var simMaxIter = sim.max_iter || 10;
        var simThreshold = sim.threshold || 'saturation';

        // ── 1. Settings Panel – collapsible, Bottom-Left ──────────────────
        var simSettingsContent = createPanel('sim_settings_panel', 'Simulation Settings',
          {bottom: '10px', left: '10px', width: '260px', borderLeft: '4px solid #3498db'});
        
        var settingsHTML =
          '<div style=\"margin-bottom:12px;\">' +
            '<label style=\"display:block; margin-bottom:5px; font-weight:bold; color:#444;\">Threshold Function</label>' +
            '<select id=\"sim_thr_sel\" style=\"width:100%; padding:4px; border-radius:4px;\">' +
              '<option value=\"saturation\"' + (simThreshold === 'saturation' ? ' selected' : '') + '>Saturation</option>' +
              '<option value=\"tanh\"' + (simThreshold === 'tanh' ? ' selected' : '') + '>Hyperbolic (Tanh)</option>' +
              '<option value=\"linear\"' + (simThreshold === 'linear' ? ' selected' : '') + '>Linear (No limit)</option>' +
            '</select>' +
          '</div>' +
          '<div style=\"margin-bottom:12px;\">' +
            '<label style=\"display:block; margin-bottom:5px; font-weight:bold; color:#444;\">Simulation Depth (Iter)</label>' +
            '<input type=\"number\" id=\"sim_depth_input\" min=\"1\" max=\"50\" value=\"' + simMaxIter + '\"' +
              ' style=\"width:100%; padding:4px; border-radius:4px; border:1px solid #ccc;\">' +
          '</div>' +
          '<div style=\"border-top:1px solid #eee; padding-top:10px;\">' +
            '<b style=\"color:#444; display:block; margin-bottom:8px;\">Scenario</b>' +
            '<div id=\"act_list\" style=\"max-height:220px; overflow-y:auto; border:1px solid #f8f9fa; padding:5px; background:#fafafa;\"></div>' +
          '</div>' +
          '<div style=\"border-top:1px solid #eee; padding-top:10px; margin-top:8px;\">' +
            '<button id=\"reset_sim\" style=\"width:100%; padding:5px 0; background:#fef9f1; border:1px solid #e0c97a; border-radius:4px; font-size:11px; font-weight:bold; color:#888; cursor:pointer;\">&#8635; Reset Scenario</button>' +
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
          bottom: '14px', left: '50%', transform: 'translateX(-50%)', width: '300px', boxSizing: 'border-box'
        });
        timelineEl.innerHTML =
          '<div style=\"display:flex; align-items:center; gap:7px; width:100%;\">' +
            '<button id=\"play_btn\" title=\"Play\" style=\"width:26px;height:26px;flex-shrink:0;border:none;border-radius:50%;background:#3498db;color:#fff;font-size:12px;cursor:pointer;padding:0;line-height:1;\">&#9654;</button>' +
            '<button id=\"pause_btn\" title=\"Pause\" style=\"width:26px;height:26px;flex-shrink:0;border:1px solid #ccc;border-radius:50%;background:#f5f5f5;color:#555;font-size:10px;cursor:pointer;padding:0;line-height:1;\">&#9646;&#9646;</button>' +
            '<input type=\"range\" id=\"sim_slider\" min=\"0\" max=\"' + simMaxIter + '\" value=\"0\"' +
              ' style=\"flex:1;accent-color:#3498db;cursor:pointer;margin:0;\">' +
            '<span id=\"iter_label\" style=\"flex-shrink:0;font-size:11px;font-weight:bold;color:#3498db;white-space:nowrap;min-width:38px;text-align:right;\">0/' + simMaxIter + '</span>' +
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
            '.wt-slider::-moz-range-track{height:6px;background:var(--wt-track-bg,#e0e0e0);border-radius:3px;}';
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
          window._simHistory = h;
          // Sync slider range
          var sliderEl = timelineEl.querySelector('#sim_slider');
          sliderEl.max = m;
          if(window._simCurrentI > m) window._simCurrentI = m;
          refreshNodes();
        };

        var updateIteration = function(idx) {
          window._simCurrentI = idx;
          if(!window._simHistory[idx]) return;
          timelineEl.querySelector('#sim_slider').value = idx;
          timelineEl.querySelector('#iter_label').innerText = idx + '/' + (parseInt(simSettingsContent.querySelector('#sim_depth_input').value) || simMaxIter);
          refreshNodes();
        };

        // ── Event listeners ───────────────────────────────────────────────
        simSettingsContent.querySelector('#sim_thr_sel').onchange = runSimulation;
        simSettingsContent.querySelector('#sim_depth_input').onchange = runSimulation;

        timelineEl.querySelector('#sim_slider').oninput = function() { updateIteration(parseInt(this.value)); };

        var simTimer = null;
        var _playStep = function(m) {
          if(window._simCurrentI >= m) { simTimer = null; return; }
          var from = window._simCurrentI;
          var to   = from + 1;
          // Update slider + label immediately
          timelineEl.querySelector('#sim_slider').value = to;
          timelineEl.querySelector('#iter_label').innerText = to + '/' + m;
          // Animate the transition
          tweenToIteration(from, to, 450);
          window._simCurrentI = to;
          simTimer = setTimeout(function() { _playStep(m); }, 550);
        };
        timelineEl.querySelector('#play_btn').onclick = function() {
          if(simTimer) { clearTimeout(simTimer); simTimer = null; }
          if(_tweenRAF) { cancelAnimationFrame(_tweenRAF); _tweenRAF = null; }
          var m = parseInt(simSettingsContent.querySelector('#sim_depth_input').value) || simMaxIter;
          if(window._simCurrentI >= m) { window._simCurrentI = 0; refreshNodes(); }
          _playStep(m);
        };
        timelineEl.querySelector('#pause_btn').onclick = function() {
          if(simTimer) { clearTimeout(simTimer); simTimer = null; }
          if(_tweenRAF) { cancelAnimationFrame(_tweenRAF); _tweenRAF = null; }
        };
        simSettingsContent.querySelector('#reset_sim').onclick = function() {
          if(simTimer) clearTimeout(simTimer);
          targetSelf = [...sim.initial_self];
          actList.querySelectorAll('.target-slider').forEach(function(s, idx) {
            var initV  = parseFloat(s.getAttribute('data-init'));
            var idealV = parseFloat(s.getAttribute('data-ideal'));
            s.value = targetSelf[idx];
            s.parentNode.querySelector('.val-badge').innerText = targetSelf[idx].toFixed(2);
            updateSliderTrack(s, initV, idealV); // reset fill color to grey
          });
          runSimulation();
          updateIteration(0);
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

        runSimulation();
      }


      // --- Event Listeners ---
      visContent.querySelector('#btn_reset').onclick = function() {
        visContent.querySelector('#dist_slider').value = 1.0;
        visContent.querySelector('#size_slider').value = 1.0;
        visContent.querySelector('#text_size_slider').value = 20;
        currentDistMult = 1.0;
        currentSizeMult = 1.0;
        currentTextSize = 20;
        refreshNodes();
        network.stabilize();
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

      visContent.querySelector('#dist_slider').oninput = function() {
        var newMult = parseFloat(this.value);
        var ratio = newMult / currentDistMult;
        currentDistMult = newMult;
        var nodesDS = network.body.data.nodes;
        var nodes = nodesDS.get();
        var updates = nodes.map(n => ({id: n.id, x: n.x * ratio, y: n.y * ratio}));
        nodesDS.update(updates);
        refreshNodes();
        network.fit({animation: false});
      };

      visContent.querySelector('#size_slider').oninput = function() {
        currentSizeMult = parseFloat(this.value);
        refreshNodes();
      };

      visContent.querySelector('#text_size_slider').oninput = function() {
        currentTextSize = parseInt(this.value);
        refreshNodes();
      };

      visContent.querySelector('#weight_slider').oninput = function() {
        var threshold = parseFloat(this.value);
        visContent.querySelector('#weight_val_txt').innerText = threshold.toFixed(2);
        var edges = network.body.data.edges;
        var updates = edges.get().map(edge => ({id: edge.id, hidden: Math.abs(edge.weight) < threshold}));
        edges.update(updates);
      };

      // Initial refresh to ensure centering on load
      refreshNodes();

      var nodeListDiv = visContent.querySelector('#node_list');
      var nodesDS = network.body.data.nodes;
      nodesDS.get().forEach(node => {
        var div = document.createElement('div');
        div.style.marginBottom = '4px';
        div.innerHTML = '<label style=\"cursor:pointer; display:flex; align-items:center; font-size:11px;\">' +
                        '<input type=\"checkbox\" class=\"node-check\" data-id=\"' + node.id + '\" ' + (node.hidden ? '' : 'checked') + ' style=\"margin-right:6px;\">' +
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

      var forceFit = function() { network.setSize('100%', el.style.height); network.redraw(); network.fit(); };
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
    g$x$color_palette_js <- list(
      "red/green" = c("#F52722", "#A5D610", "#999999", "#FFFF00"),
      "grey scale" = c("#808080", "#ffffff", "#f2f2f2", "#e5e5e5"),
      "colorblind" = c("#D55E00", "#0173B2", "#CC79A7", "#F0E442"),
      "pastel" = c("#f1677c", "#98FB98", "#F0F8FF", "#fcf087"),
      "dark" = c("#8B0000", "#006400", "#696969", "#DAA520"),
      "viridis" = c("#440154", "#35b779", "#31688e", "#fde725")
    )
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
