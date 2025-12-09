### HIDE FUNCTIONS ###


#' @importFrom dplyr select arrange case_when
#' @importFrom magrittr %>%
#' @importFrom igraph degree layout_with_graphopt layout_in_circle
#' @importFrom igraph layout_as_tree layout_with_mds layout_on_grid
#' @importFrom igraph betweenness closeness
#' @importFrom jsonlite toJSON


# Construct Names Helper Function ----------------------------------------------

.construct_names <- function(wimp) {
  paste(wimp$vertices$left_pole, "-", wimp$vertices$right_pole, sep = " ")
}

# Threshold Function -----------------------------------------------------------
.thr <- function(x, method = "saturation") {

  if (method == "none") {
    result <- x
  }

  if (method == "saturation") {
    if (x <= -1) {
      result <- -1
    }
    if (-1 < x && x < 1) {
      result <- x
    }
    if (x >= 1) {
      result <- 1
    }
  }

  if (method == "tanh") {
    result <- tanh(x)
  }

  return(result)
}


# Color function ---------------------------------------------------------------
.color_palette <- function(x) {

  # Valid color modes
  valid_modes <- c("red/green", "grey scale", "colorblind",
                   "pastel", "dark", "viridis")

  # Check if input is valid
  if (!x %in% valid_modes) {
    stop("Invalid color mode. Valid options are: ",
         paste(valid_modes, collapse = ", "))
  }

  # Order: discrepant, congruent, undefined , dilemmatic
  if (x == "red/green") {
    colors <- c("#F52722", "#A5D610", "grey", "yellow")
  }

  if (x == "grey scale") {
    colors <- c("#808080", "#ffffff", "#f2f2f2", "#e5e5e5")
  }

  if (x == "colorblind") {
    colors <- c("#D55E00", "#0173B2", "#CC79A7", "#F0E442")
  }

  if (x == "pastel") {
    colors <- c("#f1677c", "#98FB98", "#F0F8FF", "#fcf087")
  }

  if (x == "dark") {
    colors <- c("#8B0000", "#006400", "#696969", "#DAA520")
  }

  if (x == "viridis") {
    colors <- c("#440154", "#35b779", "#31688e", "#fde725")
  }

  colors
}

# Construct colors -------------------------------------------------------------
.construct_colors <- function(wimp, mode) {
  # Input validation
  if (is.null(wimp$vertices) || nrow(wimp$vertices) == 0) {
    stop("wimp$vertices must contain data")
  }

  col_sel <- .color_palette(mode)
  n <- nrow(wimp$vertices)
  labels <- .construct_names(wimp)

  # Extract values once
  self <- wimp$vertices$self
  ideal <- wimp$vertices$ideal

  # Handle missing values
  if (any(is.na(self)) || any(is.na(ideal))) {
    warning("Missing values found in self or ideal columns")
  }

  # Create indices
  idx_dilem <- ideal == 0 & !is.na(ideal)
  idx_undef <- self == 0 & !idx_dilem & !is.na(self)
  idx_congr <- sign(self) == sign(ideal) & self != 0 & ideal != 0 &
    !is.na(self) & !is.na(ideal)
  idx_discr <- !idx_dilem & !idx_undef & !idx_congr &
    !is.na(self) & !is.na(ideal)

  # Initialize colors and assign based on classification
  colors <- character(n)
  colors[idx_discr] <- col_sel[1]  # discrepant
  colors[idx_congr] <- col_sel[2]  # congruent
  colors[idx_undef] <- col_sel[3]  # undefined
  colors[idx_dilem] <- col_sel[4]  # dilemmatic

  # Handle any remaining NA cases
  colors[is.na(self) | is.na(ideal)] <- "grey50"

  # Return as matrix with proper row/column names
  result <- matrix(colors, ncol = 1, dimnames = list(labels, "color"))

  result
}

# Align wimp function ----------------------------------------------------------

.align_wimp <- function(wimp, exclude_dilemmatics = TRUE) {

  # Extract needed data
  ideal <- wimp$vertices$ideal
  swap_indices <- which(ideal < 0)

  # Identify dilemmatic constructs (ideal = 0) if needed
  if (exclude_dilemmatics) {
    dilemmatic_indices <- which(ideal == 0)
  } else {
    dilemmatic_indices <- integer(0)
  }

  # If no swapping needed and no dilemmatic exclusion, return unchanged
  if (length(swap_indices) == 0 && length(dilemmatic_indices) == 0) {
    return(wimp)
  }

  # Create a copy to modify
  aligned_wimp <- wimp

  # 1. VERTICES TRANSFORMATION
  if (length(swap_indices) > 0) {
    # Swap poles for constructs with negative ideal
    old_left <- aligned_wimp$vertices$left_pole[swap_indices]
    old_right <- aligned_wimp$vertices$right_pole[swap_indices]

    aligned_wimp$vertices$left_pole[swap_indices] <- old_right
    aligned_wimp$vertices$right_pole[swap_indices] <- old_left

    # Flip sign of self values for swapped constructs
    aligned_wimp$vertices$self[swap_indices] <-
      -aligned_wimp$vertices$self[swap_indices]

    # Take absolute value of ideal for swapped constructs
    aligned_wimp$vertices$ideal[swap_indices] <-
      abs(aligned_wimp$vertices$ideal[swap_indices])

  }

  # 2. GLOBAL MATRICES TRANSFORMATION
  if (length(swap_indices) > 0 && !is.null(aligned_wimp$global)) {
    # Transform weight_matrix: flip signs in rows and columns
    if (!is.null(aligned_wimp$global$weight_matrix)) {
      for (idx in swap_indices) {
        aligned_wimp$global$weight_matrix[idx, ] <-
          -aligned_wimp$global$weight_matrix[idx, ]
        aligned_wimp$global$weight_matrix[, idx] <-
          -aligned_wimp$global$weight_matrix[, idx]
      }
    }

    # Transform hypo_matrix: flip signs in rows for swapped constructs
    if (!is.null(aligned_wimp$global$hypo_matrix)) {
      for (idx in swap_indices) {
        aligned_wimp$global$hypo_matrix[idx, ] <-
          -aligned_wimp$global$hypo_matrix[idx, ]
      }

      # Update row/column names to reflect new pole orientations
      if (!is.null(rownames(aligned_wimp$global$hypo_matrix))) {
        for (i in swap_indices) {
          new_name <- paste("Totally", aligned_wimp$vertices$right_pole[i])
          if (i <= nrow(aligned_wimp$global$hypo_matrix)) {
            rownames(aligned_wimp$global$hypo_matrix)[i] <- new_name
          }
          if (i <= ncol(aligned_wimp$global$hypo_matrix)) {
            colnames(aligned_wimp$global$hypo_matrix)[i] <- new_name
          }
        }
      }
    }
  }

  # 3. EXCLUSION OF DILEMMATIC CONSTRUCTS
  if (exclude_dilemmatics && length(dilemmatic_indices) > 0) {
    # Remove dilemmatic constructs from all components
    keep_indices <- setdiff(seq_len(nrow(aligned_wimp$vertices)),
                            dilemmatic_indices)

    # Filter vertices
    aligned_wimp$vertices <- aligned_wimp$vertices[keep_indices, ]
    aligned_wimp$vertices$id <- seq_len(nrow(aligned_wimp$vertices))

    # Filter global matrices
    if (!is.null(aligned_wimp$global)) {
      if (!is.null(aligned_wimp$global$weight_matrix)) {
        aligned_wimp$global$weight_matrix <-
          aligned_wimp$global$weight_matrix[keep_indices, keep_indices]
      }
      if (!is.null(aligned_wimp$global$hypo_matrix)) {
        aligned_wimp$global$hypo_matrix <-
          aligned_wimp$global$hypo_matrix[keep_indices, keep_indices]
      }

      # Update n_constructs
      aligned_wimp$global$n_constructs <- length(keep_indices)
    }

    # Filter and update edges if present
    if (!is.null(aligned_wimp$edges)) {
      edges_to_keep <- !(aligned_wimp$edges$from %in% dilemmatic_indices |
                           aligned_wimp$edges$to %in% dilemmatic_indices)
      aligned_wimp$edges <- aligned_wimp$edges[edges_to_keep, ]

      # Remap vertex indices in edges
      id_mapping <- setNames(seq_along(keep_indices), keep_indices)
      aligned_wimp$edges$from <-
        id_mapping[as.character(aligned_wimp$edges$from)]
      aligned_wimp$edges$to <-
        id_mapping[as.character(aligned_wimp$edges$to)]

      # Update edge IDs
      aligned_wimp$edges$id <- paste0(aligned_wimp$edges$from, "t",
                                      aligned_wimp$edges$to)
    }
  }

  return(aligned_wimp)
}

# Hypothetical Situations Vector calculation -----------------------------------

.calc_hypo <- function(self, ideal) {

  if (is.na(self) || is.na(ideal)) return(NA_real_)

  if (self != 0) {
    result <- - sign(self)

  } else if (self == 0 && !(0 %in% ideal)) {
    result <- sign(ideal)

  } else if (self == 0 && (0 %in% ideal)) {
    result <- 1
  }

  result
}

# Dilemmatics detection ---------------------------------------------------
.which_dilemmatics <- function(wimp) {
  ideal <- wimp$vertices$ideal
  dil_indices <- which(ideal == 0)
  dil_indices
}


# PCSD Y-Axis label -------------------------------------------------------
.label_y <- function(infer) {
  if (infer == "self dynamics") "SELF DIFFERENTIAL"
  else if (infer == "impact dynamics") "IMPACT"
}

# Self Construct detection---------------------------
.self_poles <- function(self, left_pole, right_pole) {

  construct <- paste(left_pole, "-", right_pole)

  # Handle NA or missing values
  if (is.na(self)) return(construct)

  if (self < 0) return(left_pole)

  if (self > 0) return(right_pole)

  if (self == 0) return(construct)

  construct
}

# Plot Optimization Helper -----------------------------------------------------

.plot_optimization <- function(plot,
                               display_mode_bar = FALSE,
                               responsive = TRUE,
                               static_plot = FALSE) {

  config_options <- list(
    display_mode_bar = display_mode_bar,
    responsive = responsive,
    static_plot = static_plot,
    displaylogo = FALSE,
    mode_bar_buttons_to_remove = list("pan2d", "lasso2d", "select2d",
                                      "autoScale2d")
  )

  plot %>% plotly::config(
    displayModeBar = config_options$display_mode_bar,
    responsive = config_options$responsive,
    staticPlot = config_options$static_plot,
    displaylogo = config_options$displaylogo,
    modeBarButtonsToRemove = config_options$mode_bar_buttons_to_remove
  )
}


# Merge two wimps --------------------------------------------------------------
.merge_wimp <- function(wimp1, wimp2) {
  df1 <- data.frame(
    Construct = .construct_names(wimp1),
    lpoles = wimp1$vertices$left_pole,
    rpoles = wimp1$vertices$right_pole,
    index1 = seq_len(nrow(wimp1$vertices))
  )
  df2 <- data.frame(
    Construct = .construct_names(wimp2),
    lpoles = wimp2$vertices$left_pole,
    rpoles = wimp2$vertices$right_pole,
    index2 = seq_len(nrow(wimp2$vertices))
  )
  merge(df1, df2, by = 1:3, sort = FALSE)
}

# Compatibility merge wimps ----------------------------------------------------
.compatibility_merge_wimp <- function(wimp1, wimp2) {
  m <- nrow(.merge_wimp(wimp1, wimp2))
  n1 <- nrow(wimp1$vertices)
  n2 <- nrow(wimp2$vertices)
  if (m == 0) return("Incompatibility")
  if (m == n1 && m == n2) return("Full Compatibility")
  "Partial Compatibitily"
}

# Tversky Similarity function---------------------------------------------------
.sim_index <- function(x, y, alpha = .5, beta = .5) {

  s_vec <- ifelse(y^2 > 0.25,
                  (-(x - y)^2) / (y^2) + 1,
                  (-(x - y)^2) / ((1 - abs(y)^2) + 1))

  vec_int <- s_vec[which(x * y > 0)]
  vec_x_minus_y <- x[which(x * y <= 0)]
  vec_y_minus_x <- y[which(x * y <= 0)]

  int_xy <- sum(abs(vec_int))
  x_minus_y <- sum(abs(vec_x_minus_y))
  y_minus_x <- sum(abs(vec_y_minus_x))

  sim_ratio <- int_xy / (int_xy + alpha * x_minus_y + beta * y_minus_x)
  return(sim_ratio)

}

# Similarity function-----------------------------------------------------------
.sim <- function(s, i) {

  result <- 1 - ((s - i)^2 / 4)
  result
}

# Impact function---------------------------------------------------------------
.impact <- function(w, s, i) {
  s_next <- .thr(s + w)
  result <- .sim(s_next, i) - .sim(s, i)
  result
}

# JavaScript Hulls Generator ---------------------------------------------------
.create_js_hulls <- function(pad_side, rounding, uniq_cats, fills, strokes) {
  # nolint start
  paste0(
    "function(ctx){
  var net=this;
  var PAD=", pad_side, ";
  var ROUND=", rounding, ";
  var LINE_W=2;
  var MAX_OFFSET_DIST=150;

  var CAT = ", jsonlite::toJSON(uniq_cats, auto_unbox=TRUE), ";
  var FILL= ", jsonlite::toJSON(as.character(fills), auto_unbox=TRUE), ";
  var STROK=", jsonlite::toJSON(as.character(strokes), auto_unbox=TRUE), ";

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
      
      // Calculate angle to detect sharp corners
      var dot=e0.x*e1.x+e0.y*e1.y; var angle=Math.acos(Math.max(-1,Math.min(1,-dot)));
      var isSharpCorner=angle<Math.PI/3; // Menos de 60 grados
      
      var pB={x:p1.x+n0.x*d,y:p1.y+n0.y*d}, pC={x:p1.x+n1.x*d,y:p1.y+n1.y*d};
      var q=lineIntersect(pB,e0,pC,e1);
      
      // Validate intersection and distance
      if(q && distance(p1,q)<=MAX_OFFSET_DIST && !isSharpCorner){
        out[i]=q;
      } else {
        // Improved fallback for problematic corners
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
    ctx.save(); ctx.globalCompositeOperation='source-over';
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
        
        // Adapt rounding radius based on angle
        var dot=(v1x/l1)*(-v2x/l2)+(v1y/l1)*(-v2y/l2);
        var angle=Math.acos(Math.max(-1,Math.min(1,dot)));
        var angleRatio=Math.max(0.3,Math.sin(angle*0.5));
        var adaptiveRound=ROUND*angleRatio;
        
        var rr=Math.min(adaptiveRound,0.4*l1,0.4*l2);
        v1x/=l1; v1y/=l1; v2x/=l2; v2y/=l2;
        var p1_in={x:p1.x-v1x*rr,y:p1.y-v1y*rr}, p1_out={x:p1.x+v2x*rr,y:p1.y+v2y*rr};
        if(i===0) ctx.moveTo(p1_in.x,p1_in.y); else ctx.lineTo(p1_in.x,p1_in.y);
        
        // Usar cuadraticCurveTo para esquinas muy agudas
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
    var xs = pts.map(p => p.x), ys = pts.map(p => p.y);
    var minX = Math.min.apply(null,xs), maxX = Math.max.apply(null,xs);
    var minY = Math.min.apply(null,ys);
    var cx = (minX + maxX) / 2, y = minY - 12;
    ctx.save(); ctx.font='bold 14px sans-serif'; ctx.textAlign='center'; ctx.textBaseline='bottom';
    ctx.strokeStyle = stroke; ctx.lineWidth = 4; ctx.strokeText(text, cx, y);
    ctx.fillStyle = 'rgba(0,0,0,0.90)'; ctx.fillText(text, cx, y); ctx.restore();
  }

  var nodes=net.body.data.nodes.get(); var bycat={};
  nodes.forEach(function(n){
    if(n.hidden) return;
    if(n.category==null) return;
    var c=String(n.category);
    if(c===''||c==='NA') return;
    (bycat[c]||(bycat[c]=[])).push(n.id);
  });

  for(var i=0;i<CAT.length;i++){
    var cat = CAT[i], ids = bycat[cat];
    if(!ids || ids.length===0) continue;
    var pos = net.getPositions(ids);
    var pts = ids.map(function(id){ return {x:pos[id].x, y:pos[id].y}; });

    var outline;
    if(pts.length<=2){
      // Mejorado: padding adaptativo para pocos puntos
      var xs=pts.map(p => p.x), ys=pts.map(p => p.y);
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
      // Mejorado: adaptive padding based on point dispersion
      var xs=pts.map(p => p.x), ys=pts.map(p => p.y);
      var avgX=xs.reduce((a,b) => a+b)/xs.length, avgY=ys.reduce((a,b) => a+b)/ys.length;
      var avgDist=pts.reduce((sum,p) => sum+distance({x:avgX,y:avgY},p),0)/pts.length;
      var adaptivePad=Math.max(PAD, avgDist*0.15+15);
      
      var hull=convexHull(pts); 
      outline=offsetConvex(hull, adaptivePad);
      
      // Fallback si offset falla
      if(!outline){
        var minX=Math.min.apply(null,xs)-adaptivePad, maxX=Math.max.apply(null,xs)+adaptivePad;
        var minY=Math.min.apply(null,ys)-adaptivePad, maxY=Math.max.apply(null,ys)+adaptivePad;
        outline=[{x:minX,y:minY},{x:maxX,y:minY},{x:maxX,y:maxY},{x:minX,y:maxY}];
      }
    }

    draw(ctx, outline, FILL[i], STROK[i]);
    drawLabelAbove(ctx, outline, String(cat), STROK[i]);
  }
}"
  )
  # nolint end
}

# Salience Parameter Estimation ------------------------------------------------

# Calculate structural coefficients for SSI (Definition 9)
.calc_structural_coefs <- function(wimp) {
  
  self <- wimp$vertices$self
  ideal <- wimp$vertices$ideal
  
  # Helper function: local similarity g(s_i, d_i) from Definition 5
  g_similarity <- function(s, d) {
    numerator <- (s - d)^2
    denominator <- pmax(d^2, (1 - abs(d))^2)
    1 - (numerator / denominator)
  }
  
  # Definition 4: Set relations
  # (i) Shared attributes: s_i · d_i > 0 (congruent)
  congruent_idx <- which(self * ideal > 0)
  
  # (ii) Discrepancies: s_i · d_i ≤ 0 (S \ I)
  discrepant_idx <- which(self * ideal <= 0)
  
  # Definition 5: Magnitude of shared attributes f(S ∩ I)
  f_shared <- if (length(congruent_idx) > 0) {
    sum(g_similarity(self[congruent_idx], ideal[congruent_idx]))
  } else {
    1e-6  # Avoid division by zero
  }
  
  # Definition 6: Magnitude of self-now discrepancies f(S \ I)
  f_discrepancy <- if (length(discrepant_idx) > 0) {
    sum(abs(self[discrepant_idx]))
  } else {
    0
  }
  
  # Definition 7: Magnitude of aspirational gaps f(I \ S)
  f_aspiration <- if (length(discrepant_idx) > 0) {
    sum(abs(ideal[discrepant_idx]))
  } else {
    0
  }
  
  # Definition 9: Structural coefficients
  omega_alpha <- f_discrepancy / f_shared
  omega_beta <- f_aspiration / f_shared
  
  list(
    omega_alpha = omega_alpha,
    omega_beta = omega_beta,
    f_shared = f_shared,
    f_discrepancy = f_discrepancy,
    f_aspiration = f_aspiration
  )
}

.estimate_ssi_parameters <- function(wimp) {
  
  # --- 1. Extracción y Limpieza ---
  self <- wimp$vertices$self
  ideal <- wimp$vertices$ideal
  preference <- wimp$vertices$preference
  
  if (is.character(preference)) preference <- as.numeric(preference)
  if (is.null(preference) || all(is.na(preference))) preference <- rep(0, length(self))
  
  # --- 2. Clasificación de Constructos ---
  # Congruencia: Signos coinciden
  is_congruent <- (sign(self) * sign(ideal)) > 0
  idx_cong <- which(is_congruent)
  idx_disc <- which(!is_congruent)
  
  # --- 3. Cálculo de Importancias (Weighted Saliency) ---
  # A) Congruencia
  w_cong_values <- numeric(length(idx_cong))
  if (length(idx_cong) > 0) {
    matches <- sign(preference[idx_cong]) == sign(self[idx_cong])
    w_cong_values[matches] <- abs(preference[idx_cong][matches])
    # Importante: Si no coincide, asumimos un valor residual pequeño en vez de 0 absoluto
    # para evitar varianzas colapsadas
    w_cong_values[!matches] <- 0.1 
  }
  
  # B) Discrepancia Yoica (Alpha)
  w_disc_values <- numeric(length(idx_disc))
  if (length(idx_disc) > 0) {
    matches_self <- sign(preference[idx_disc]) == sign(self[idx_disc])
    w_disc_values[matches_self] <- abs(preference[idx_disc][matches_self])
    w_disc_values[!matches_self] <- 0.1
  }
  
  # C) Aspiración (Beta)
  w_asp_values <- numeric(length(idx_disc))
  if (length(idx_disc) > 0) {
    matches_ideal <- sign(preference[idx_disc]) == sign(ideal[idx_disc])
    w_asp_values[matches_ideal] <- abs(preference[idx_disc][matches_ideal])
    w_asp_values[!matches_ideal] <- 0.1
  }
  
  # --- 4. Estadísticos con Suavizado ---
  # Usamos un suavizado bayesiano simple (añadir pseudo-observaciones) 
  # para que las medias nunca sean 0 y las varianzas no exploten.
  
  safe_mean <- function(x) {
    if(length(x) == 0) return(0.5)
    mean(c(x, 1), na.rm=TRUE) # "Push" suave hacia 1 para evitar ceros
  }
  
  safe_var <- function(x) {
    if(length(x) < 2) return(0.25) # Varianza por defecto alta si hay pocos datos
    var(c(x, 0.5, 1), na.rm=TRUE) # Añadir variabilidad artificial para estabilidad
  }

  w_bar_cong <- safe_mean(w_cong_values)
  w_bar_disc <- safe_mean(w_disc_values)
  w_bar_asp  <- safe_mean(w_asp_values)
  
  var_cong <- safe_var(w_cong_values)
  var_disc <- safe_var(w_disc_values)
  var_asp  <- safe_var(w_asp_values)
  
  # --- 5. Parámetros Normalizados ---
  denom_alpha <- w_bar_disc + w_bar_cong
  denom_beta  <- w_bar_asp + w_bar_cong
  
  mu_alpha <- w_bar_disc / denom_alpha
  mu_beta  <- w_bar_asp / denom_beta
  
  # --- 6. Varianza (Método Delta) ---
  deriv_disc_a <- w_bar_cong / (denom_alpha^2)
  deriv_cong_a <- -w_bar_disc / (denom_alpha^2)
  sigma2_alpha <- (deriv_disc_a^2 * var_disc) + (deriv_cong_a^2 * var_cong)
  
  deriv_asp_b  <- w_bar_cong / (denom_beta^2)
  deriv_cong_b <- -w_bar_asp / (denom_beta^2)
  sigma2_beta  <- (deriv_asp_b^2 * var_asp) + (deriv_cong_b^2 * var_cong)
  
  # --- 7. SANITY CHECK (La clave del arreglo) ---
  # Forzamos que la varianza esté en un rango razonable para visualización.
  # Una varianza > 0.1 en una escala 0-1 aplana demasiado la curva.
  # Una varianza < 0.001 hace un pico invisible.
  
  clamp <- function(x, min_val=0.005, max_val=0.05) {
    max(min(x, max_val), min_val)
  }
  
  sigma2_alpha <- clamp(sigma2_alpha)
  sigma2_beta  <- clamp(sigma2_beta)
  
  # --- 8. Función PDF ---
  pdf_function <- function(alpha_grid, beta_grid) {
    const <- 1 / (2 * pi * sqrt(sigma2_alpha * sigma2_beta))
    z_score <- ((alpha_grid - mu_alpha)^2 / sigma2_alpha) + 
               ((beta_grid - mu_beta)^2  / sigma2_beta)
    return(const * exp(-0.5 * z_score))
  }
  
  list(
    mu_alpha = mu_alpha,
    mu_beta = mu_beta,
    sigma2_alpha = sigma2_alpha,
    sigma2_beta = sigma2_beta,
    pdf_function = pdf_function
  )
}
