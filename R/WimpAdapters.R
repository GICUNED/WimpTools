## ADAPTER FUNCTIONS FOR NEW WIMP FORMAT ##
## These functions handle both old and new wimp structures transparently

# Detect if wimp is new format (has global, vertices, edges) or old format
.wimp_is_new <- function(wimp) {
  return(!is.null(wimp$global) && !is.null(wimp$vertices))
}

# Get number of constructs
.wimp_n_constructs <- function(wimp) {
  if(.wimp_is_new(wimp)) {
    if(is.data.frame(wimp$vertices)) return(nrow(wimp$vertices))
    if(is.list(wimp$vertices)) return(length(wimp$vertices$id))
  }
  stop("Invalid wimp: expected new format with $global and $vertices")
}

# Get vertices dataframe (new format) or reconstruct from old
.wimp_get_vertices <- function(wimp) {
  if(.wimp_is_new(wimp)) return(wimp$vertices)
  stop("Invalid wimp: expected new format with $vertices")
}

# Get edges dataframe (new format) or reconstruct from old
.wimp_get_edges <- function(wimp) {
  if(.wimp_is_new(wimp)) return(wimp$edges)
  stop("Invalid wimp: expected new format with $edges")
}

# Get weight matrix (n x n standardized implications)
.wimp_get_weights_matrix <- function(wimp) {
  if(is.list(wimp$global) && !is.null(wimp$global$wmatrix)) return(wimp$global$wmatrix)
  # Fallback: reconstruct from edges
  edges <- .wimp_get_edges(wimp)
  n <- as.numeric(.wimp_n_constructs(wimp))
  wmatrix <- matrix(0, nrow = n, ncol = n)
  if(!is.null(edges) && nrow(edges) > 0) {
    for(r in seq_len(nrow(edges))) {
      i <- as.integer(edges[r, "from"]); j <- as.integer(edges[r, "to"])
      wmatrix[i, j] <- as.numeric(edges[r, "weight"])
    }
  }
  return(wmatrix)
}

# Get standardized self vector
.wimp_get_self <- function(wimp) {
  return(wimp$vertices$self)
}

# Get standardized ideal vector
.wimp_get_ideal <- function(wimp) {
  return(wimp$vertices$ideal)
}

# Get left poles vector
.wimp_get_left_poles <- function(wimp) {
  return(wimp$vertices$lpole)
}

# Get right poles vector
.wimp_get_right_poles <- function(wimp) {
  return(wimp$vertices$rpole)
}

# Get construct names (paste lpole - rpole)
.wimp_get_construct_names <- function(wimp) {
  lpoles <- wimp$vertices$lpole
  rpoles <- wimp$vertices$rpole
  return(paste(lpoles, "-", rpoles, sep = " "))
}

# Get scale min/max
.wimp_get_scale <- function(wimp) {
  if(is.list(wimp$global) && !is.null(wimp$global$scale)){
    smin <- as.numeric(wimp$global$scale[1])
    smax <- as.numeric(wimp$global$scale[2])
    return(list(min = smin, max = smax, center = (smin + smax) / 2))
  }
  stop("Scale information not available in wimp$global")
}

# Get implications matrix (before weight calculation)
.wimp_get_implications_matrix <- function(wimp) {
  stop("Implications matrix not available in new format. Use weight matrix instead.")
}

# Get hypothetical matrix for PCSD-type analyses
.wimp_get_hypothetical_matrix <- function(wimp) {
  return(.wimp_get_weights_matrix(wimp))
}

# Reconstruct minimal legacy structure (for backward compatibility)
.wimp_as_legacy <- function(wimp) {
  stop("Legacy conversion not supported; use new wimp format only.")
}

# Get construct indices (congruents, discrepants, dilemmatics, undefined)
.wimp_get_construct_indices <- function(wimp) {
  self_vec <- .wimp_get_self(wimp)
  ideal_vec <- .wimp_get_ideal(wimp)
  
  congruents <- which(sign(self_vec) == sign(ideal_vec) & self_vec != 0 & ideal_vec != 0)
  discrepants <- which(sign(self_vec) != sign(ideal_vec) & self_vec != 0 & ideal_vec != 0)
  dilemmatics <- which(ideal_vec == 0)
  undefined <- setdiff(seq_along(self_vec), c(congruents, discrepants, dilemmatics))
  
  return(list(
    congruents = congruents,
    discrepants = discrepants,
    dilemmatics = dilemmatics,
    undefined = undefined
  ))
}
