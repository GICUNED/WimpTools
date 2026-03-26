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
#' @param height Character string specifying the graph height. Default is
#'        "700px".
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
                    height = "700px", color = "red/green", layout = "graphopt",
                    show = TRUE, hide_direct = FALSE,
                    areas = FALSE, area_attr = "category", area_color = NA,
                    pad_side = 50, rounding = 10, min_weight = 0,
                    interactive_options = TRUE) {

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

  # Create vertex properties
  vertex_name <- .create_vertex_names(vertex_vector, lpoles, rpoles)
  congruency <- .calculate_congruency(vertex_vector, ideal_vector, color)
  vertex_vadjust <- sapply(abs(vertex_vector),
                           function(x) -4 * x^2 - 25 * x - 42)

  # Build vertices data frame
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
    color = congruency$color,
    orig_color = congruency$color,
    self_val = as.numeric(vertex_vector),
    ideal_val = as.numeric(ideal_vector),
    shadow = TRUE,
    hidden = !show,
    font.size = 20,
    font.strokeWidth = 3,
    font.vadjust = vertex_vadjust
  )
  vertex[[area_attr]] <- as.character(area_vec)

  # Internal helper for layout pre-calculation
  .get_all_layouts <- function(wmatrix, vertex, area_attr) {
    ig <- igraph::graph_from_adjacency_matrix(wmatrix, weight = TRUE, mode = "directed")
    
    # igraph layouts
    # Normalize coordinates to a consistent range [-500, 500]
    .scale <- function(m) {
      if (nrow(m) == 0) return(m)
      # Normalize to [0, 1] then to [-500, 500]
      for (i in 1:2) {
        rng <- range(m[, i])
        span <- rng[2] - rng[1]
        if (span < 1e-9) span <- 1
        m[, i] <- (m[, i] - rng[1]) / span - 0.5
      }
      m * 850
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
    
    # Add "Original" layout (pre-calculate based on parameter)
    # This ensures we can always return to the initial state
    init_layout_name <- .convert_layout_name(layout)
    if (init_layout_name == "layout_as_tree") {
       layouts[["original"]] <- .scale(igraph::layout_as_tree(ig, circular = TRUE))
    } else if (init_layout_name == "layout_in_circle") {
       layouts[["original"]] <- .scale(igraph::layout_in_circle(ig))
    } else if (init_layout_name == "layout_with_mds") {
       layouts[["original"]] <- .scale(igraph::layout_with_mds(ig))
    } else if (init_layout_name == "layout_on_grid") {
       layouts[["original"]] <- .scale(igraph::layout_on_grid(ig))
    } else if (init_layout_name == "areas") {
       # already handled if areas=TRUE
    } else {
       layouts[["original"]] <- .scale(igraph::layout_with_graphopt(ig))
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
    
    edges <- data.frame(
      from = as.character(edges_raw$from),
      to = as.character(edges_raw$to),
      width = 2 * abs(edges_raw$weight),
      arrows = "to",
      dashes = edge_props$dashes,
      smooth = edge_curved,
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
      "grid"     = "layout_on_grid"
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
      visInteraction(navigationButtons = TRUE, multiselect = TRUE) %>%
      visPhysics(enabled = FALSE)
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
      visInteraction(navigationButtons = TRUE, multiselect = TRUE)
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
      visInteraction(navigationButtons = TRUE, multiselect = TRUE)
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
      var network = this.network;
      var panel = document.createElement('div');
      panel.className = 'wimp-options-panel';
      panel.style.position = 'absolute';
      panel.style.top = '10px';
      panel.style.right = '10px';
      panel.style.zIndex = '1000';
      panel.style.backgroundColor = 'rgba(255, 255, 255, 0.95)';
      panel.style.padding = '10px';
      panel.style.borderRadius = '8px';
      panel.style.boxShadow = '0 2px 15px rgba(0,0,0,0.15)';
      panel.style.border = '1px solid #ddd';
      panel.style.fontFamily = 'Segoe UI, Tahoma, Geneva, Verdana, sans-serif';
      panel.style.fontSize = '12px';
      panel.style.maxHeight = '90%';
      panel.style.overflowY = 'auto';
      panel.style.width = '100px';
      panel.style.transition = 'all 0.3s ease';
      
      el.appendChild(panel);

      // --- Header and Minimize Logic ---
      var header = document.createElement('div');
      header.style.display = 'flex';
      header.style.justifyContent = 'space-between';
      header.style.alignItems = 'center';
      header.style.cursor = 'pointer';
      header.style.userSelect = 'none';
      
      header.innerHTML = '<span style=\"font-weight:bold; color:#2c3e50;\">Options</span>' +
                         '<span id=\"toggle_btn\" style=\"font-weight:bold; color:#95a5a6; font-size:16px; width:20px; textAlign:center;\">+</span>';
      panel.appendChild(header);

      var content = document.createElement('div');
      content.id = 'panel_content';
      content.style.display = 'none'; 
      content.style.marginTop = '10px';
      content.style.borderTop = '1px solid #eee';
      content.style.paddingTop = '10px';
      panel.appendChild(content);

      header.onclick = function() {
        var isHidden = content.style.display === 'none';
        content.style.display = isHidden ? 'block' : 'none';
        header.querySelector('#toggle_btn').innerText = isHidden ? '−' : '+';
        panel.style.width = isHidden ? '220px' : '100px';
        panel.style.padding = isHidden ? '12px' : '10px';
      };      // --- Palette Section ---
      var paletteDiv = document.createElement('div');
      paletteDiv.style.marginBottom = '15px';
      
      var pOpts = ['red/green', 'grey scale', 'colorblind', 'pastel', 'dark', 'viridis'];
      var pLabels = ['Red / Green', 'Grey Scale', 'Colorblind', 'Pastel', 'Dark', 'Viridis'];
      var pSelectHtml = '<div style=\"margin-bottom:5px; font-weight:bold; color:#444;\">Color Palette</div>' +
                         '<select id=\"palette_sel\" style=\"width:100%; padding:4px; font-size:11px;\">';
      
      for(var i=0; i<pOpts.length; i++){
        var sel = (pOpts[i] === x.initial_palette) ? ' selected=\"selected\"' : '';
        pSelectHtml += '<option value=\"' + pOpts[i] + '\"' + sel + '>' + pLabels[i] + '</option>';
      }
      pSelectHtml += '</select>';
      
      paletteDiv.innerHTML = pSelectHtml;
      content.appendChild(paletteDiv);
      
      var getPaletteColor = function(v, i, scheme) {
        var x_val = v / i;
        var colors = {
          'red/green':  [\"#F52722\", \"#A5D610\", \"#999999\", \"#FFFF00\"],
          'grey scale': [\"#808080\", \"#ffffff\", \"#f2f2f2\", \"#e5e5e5\"],
          'colorblind': [\"#D55E00\", \"#0173B2\", \"#CC79A7\", \"#F0E442\"],
          'pastel':     [\"#f1677c\", \"#98FB98\", \"#F0F8FF\", \"#fcf087\"],
          'dark':       [\"#8B0000\", \"#006400\", \"#696969\", \"#DAA520\"],
          'viridis':    [\"#440154\", \"#35b779\", \"#31688e\", \"#fde725\"]
        };
        var p = colors[scheme] || colors['red/green'];
        if (x_val < 0 && x_val !== -Infinity) return p[0];
        if (x_val > 0 && x_val !== Infinity) return p[1];
        if (x_val === 0) return p[2];
        return p[3];
      };
 
      var updatePalette = function(scheme) {
        var nodesDS = network.body.data.nodes;
        var nodesUpdates = nodesDS.get().map(function(node) {
          var c;
          if (scheme === x.initial_palette) {
            c = node.orig_color;
          } else {
            c = getPaletteColor(node.self_val, node.ideal_val, scheme);
          }
          return {id: node.id, color: c};
        });
        nodesDS.update(nodesUpdates);
        
        var edgesDS = network.body.data.edges;
        var edgesUpdates = edgesDS.get().map(function(edge) {
          if(scheme === 'grey scale') {
             return {id: edge.id, color: '#999999', dashes: edge.weight < 0};
          } else {
             return {id: edge.id, color: edge.orig_color, dashes: edge.orig_dashes}; 
          }
        });
        edgesDS.update(edgesUpdates);
      };
      
      // Initialize with current selection
      updatePalette(paletteDiv.querySelector('#palette_sel').value);

      paletteDiv.querySelector('#palette_sel').onchange = function() {
        updatePalette(this.value);
      };      // --- Layout Section ---
      var layoutDiv = document.createElement('div');
      layoutDiv.style.marginBottom = '15px';
      var layoutOptions = '<option value=\"original\">Original</option>' +
                          '<option value=\"graphopt\">GraphOpt</option>' +
                          '<option value=\"circle\">Circle</option>' +
                          '<option value=\"mds\">MDS</option>' +
                          '<option value=\"grid\">Grid</option>' +
                          '<option value=\"tree\">Tree (Circular)</option>';
      
      if (x.layouts.areas) {
        layoutOptions += '<option value=\"areas\">Areas</option>';
      }

      layoutDiv.innerHTML = '<div style=\"margin-bottom:5px; font-weight:bold; color:#444;\">Layout</div>' +
                            '<select id=\"layout_sel\" style=\"width:100%; padding:4px; font-size:11px; margin-bottom:5px;\">' +
                            layoutOptions +
                            '</select>';
      content.appendChild(layoutDiv);
      
      layoutDiv.querySelector('#layout_sel').onchange = function() {
        var layoutKey = this.value;
        var layoutData = x.layouts[layoutKey];
        if (layoutData) {
          network.setOptions({
            physics: {enabled: false},
            edges: {smooth: {type: 'curvedCW', roundness: 0.15}}
          });
          var nodesDS = network.body.data.nodes;
          var updates = [];
          for (var i = 0; i < layoutData.id.length; i++) {
            updates.push({
              id: String(layoutData.id[i]), 
              x: layoutData.x[i], 
              y: layoutData.y[i]
            });
          }
          nodesDS.update(updates);
          network.fit({animation: true});
        }
      };

      // --- Weight Filter Section ---
      var sliderSection = document.createElement('div');
      sliderSection.style.marginBottom = '15px';
      sliderSection.style.borderTop = '1px solid #eee';
      sliderSection.style.paddingTop = '10px';
      
      var maxW = x.max_weight.toFixed(2);
      sliderSection.innerHTML = '<div style=\"margin-bottom:8px; font-weight:bold; color:#444;\">Edge Weight Filter</div>' +
                                '<input type=\"range\" id=\"min_weight_slider\" min=\"0\" max=\"' + x.max_weight + '\" step=\"0.01\" value=\"' + (x.min_weight || 0) + '\" style=\"width:100%; cursor:pointer;\">' +
                                '<div style=\"margin-top:6px; display:flex; justify-content:space-between; font-family:monospace;\">' +
                                '<span>Min: <b id=\"weight_val\" style=\"color:#2c3e50;\">' + (x.min_weight || 0).toFixed(2) + '</b></span>' +
                                '<span style=\"color:#7f8c8d;\">M: ' + maxW + '</span></div>';
      content.appendChild(sliderSection);
      
      var slider = sliderSection.querySelector('#min_weight_slider');
      slider.addEventListener('input', function() {
        var threshold = parseFloat(this.value);
        sliderSection.querySelector('#weight_val').innerText = threshold.toFixed(2);
        var edges = network.body.data.edges;
        var updates = edges.get().map(function(edge) {
          return {id: edge.id, hidden: Math.abs(edge.weight) < threshold};
        });
        edges.update(updates);
      });

      // --- Node Filter Section ---
      var nodeSection = document.createElement('div');
      nodeSection.innerHTML = '<div style=\"margin-bottom:8px; font-weight:bold; color:#444;\">Constructs Checklist</div>' +
                               '<div style=\"margin-bottom:8px; display:flex; gap:5px;\">' +
                               '<button id=\"check_all\" style=\"flex:1; cursor:pointer; font-size:10px; padding:2px;\">All</button>' +
                               '<button id=\"uncheck_all\" style=\"flex:1; cursor:pointer; font-size:10px; padding:2px;\">None</button>' +
                               '</div>';
      
      var list = document.createElement('div');
      list.style.maxHeight = '200px';
      list.style.overflowY = 'auto';
      list.style.paddingRight = '5px';
      list.style.borderTop = '1px solid #f0f0f0';
      list.style.paddingTop = '8px';
      
      var nodesDS = network.body.data.nodes;
      nodesDS.get().forEach(function(node) {
         var item = document.createElement('div');
         item.style.marginBottom = '6px';
         var checked = node.hidden ? '' : 'checked';
         item.innerHTML = '<label style=\"cursor:pointer; display:flex; align-items:flex-start; line-height:1.2; font-size:11px;\">' +
                          '<input type=\"checkbox\" class=\"node-check\" data-id=\"' + node.id + '\" ' + checked + ' style=\"margin-top:1px; margin-right:8px;\">' +
                          '<span style=\"word-break: break-word;\">' + (node.label || 'Node ' + node.id) + '</span>' +
                          '</label>';
         list.appendChild(item);
      });
      
      nodeSection.appendChild(list);
      content.appendChild(nodeSection);
      
      var updateNodeVisibility = function(id, isVisible) {
        var nodeObj = nodesDS.get(id);
        if (nodeObj) {
          nodeObj.hidden = !isVisible;
          nodesDS.update(nodeObj);
        }
      };

      list.addEventListener('change', function(e) {
        if (e.target.classList.contains('node-check')) {
          var nodeId = String(e.target.getAttribute('data-id'));
          updateNodeVisibility(nodeId, e.target.checked);
        }
      });
      
      content.querySelector('#check_all').onclick = function() {
         list.querySelectorAll('.node-check').forEach(function(c) { 
           var nodeId = String(c.getAttribute('data-id'));
           c.checked = true; 
           updateNodeVisibility(nodeId, true);
         });
      };

      content.querySelector('#uncheck_all').onclick = function() {
         list.querySelectorAll('.node-check').forEach(function(c) { 
           var nodeId = String(c.getAttribute('data-id'));
           c.checked = false; 
           updateNodeVisibility(nodeId, false);
         });
      };
    }
    "
    g$x$interactive_options <- interactive_options
    g$x$min_weight    <- min_weight
    g$x$max_weight    <- max_w
    g$x$layouts       <- .get_all_layouts(wmatrix, vertex, area_attr)
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

  # Select the vertex vector for the requested iteration
  vertex_vector <- scn[[1]][niter + 1, ]
  digraph(wimp = scn, vertex_vector = vertex_vector, ...)
}
