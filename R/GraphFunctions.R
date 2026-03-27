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
#' @param width Character string specifying the graph width. Default is "100%".
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
#'   \item \strong{Node sizes}: Proportional to absolute self values
#'   \item \strong{Edge colors}: Represent relationship types
#'   \item \strong{Edge directions}: Show influence flow between constructs
#' }
#'
#' @author Alejandro Sanfeliciano
#'
#' @return A \code{visNetwork} interactive graph object.
#'
#' @import visNetwork
#' @importFrom htmlwidgets JS onRender
#' @importFrom jsonlite toJSON
#' @importFrom magrittr %>%
#' @importFrom visNetwork visNetwork visOptions visInteraction visPhysics
#'             visIgraphLayout visEvents
#' @export
digraph <- function(wimp, vertex_vector = NA, ideal_vector = NA, width = "100%",
                    height = "90vh", color = "red/green", layout = "graphopt",
                    show = TRUE, hide_direct = FALSE,
                    areas = FALSE, area_attr = "category", area_color = NA,
                    pad_side = 50, rounding = 10, min_weight = 0,
                    interactive_options = TRUE, sim_data = NULL, ...) {

  # ==========================================
  # INPUT VALIDATION
  # ==========================================
  allowed_colors  <- c("red/green", "grey scale", "colorblind", "pastel", "dark", "viridis")
  allowed_layouts <- c("graphopt", "circle", "rtcircle", "tree", "mds", "grid", "areas")

  if (!(inherits(wimp, "wimp") || inherits(wimp, "scn"))) stop("'wimp' must be 'wimp' or 'scn' object.")

  n_constructs <- if (inherits(wimp, "wimp")) nrow(wimp$vertices) else length(wimp$constructs[[1]])

  if (!is.na(vertex_vector[1])) {
    if (length(vertex_vector) != n_constructs) stop("Length mismatch for vertex_vector.")
  }
  if (!is.na(ideal_vector[1])) {
    if (length(ideal_vector) != n_constructs) stop("Length mismatch for ideal_vector.")
  }

  # ==========================================
  # HELPER FUNCTIONS
  # ==========================================

  .extract_wimp_data <- function(wimp, area_attr) {
    if (inherits(wimp, "wimp")) {
      wmatrix <- if (!is.null(wimp$global$weight_matrix)) wimp$global$weight_matrix else matrix(0, nrow(wimp$vertices), nrow(wimp$vertices))
      area_vec <- if (area_attr %in% names(wimp$vertices)) wimp$vertices[[area_attr]] else rep(NA_character_, nrow(wimp$vertices))
      return(list(lpoles = wimp$vertices$left_pole, rpoles = wimp$vertices$right_pole, wmatrix = wmatrix, self = wimp$vertices$self, ideal = wimp$vertices$ideal, area_vec = area_vec))
    } else if (inherits(wimp, "scn")) {
      return(list(lpoles = wimp$constructs[[1]], rpoles = wimp$constructs[[2]], wmatrix = wimp$weights, self = wimp[["self"]][[1]], ideal = wimp[["self"]][[2]], area_vec = rep(NA_character_, length(wimp$constructs[[1]]))))
    }
  }

  .calculate_congruency <- function(vertex_vector, ideal_vector, color) {
    congruency_vector <- vertex_vector / ideal_vector
    vertex_color <- sapply(congruency_vector, function(x) {
      if (is.na(x) || is.nan(x)) return(.color_palette(color)[3])
      if (is.infinite(x)) return(.color_palette(color)[4])
      if (x < 0) return(.color_palette(color)[1])
      if (x > 0) return(.color_palette(color)[2])
      return(.color_palette(color)[3])
    })
    vertex_group <- sapply(congruency_vector, function(x) {
      if (is.na(x) || is.nan(x)) return("Undefined")
      if (is.infinite(x)) return("Dilemmatic")
      if (x < 0) return("Discrepant")
      if (x > 0) return("Congruent")
      return("Undefined")
    })
    list(color = vertex_color, group = vertex_group)
  }

  .extract_edges <- function(wmatrix) {
    idx <- which(wmatrix != 0, arr.ind = TRUE)
    if (length(idx) == 0) return(data.frame(from=character(0), to=character(0), weight=numeric(0)))
    data.frame(from = idx[,1], to = idx[,2], weight = wmatrix[idx])
  }

  # ==========================================
  # MAIN LOGIC
  # ==========================================
  wd <- .extract_wimp_data(wimp, area_attr)
  if (is.na(vertex_vector[1])) vertex_vector <- wd$self
  if (is.na(ideal_vector[1])) ideal_vector <- wd$ideal
  vertex_vector <- sapply(vertex_vector, .thr)

  signs <- sign(vertex_vector); signs[signs == 0] <- 1
  wmatrix <- diag(signs, nrow = length(signs)) %*% wd$wmatrix %*% diag(signs, nrow = length(signs))

  poles <- paste(wd$lpoles, "-", wd$rpoles)
  vnames <- sapply(seq_along(vertex_vector), function(i) if(vertex_vector[i]<0) wd$lpoles[i] else if(vertex_vector[i]>0) wd$rpoles[i] else poles[i])
  cg <- .calculate_congruency(vertex_vector, ideal_vector, color)

  vertex <- data.frame(
    id = as.character(seq_along(vertex_vector)), label = vnames, group = cg$group,
    size = 30 * abs(vertex_vector) + 20, shape = "dot",
    color = cg$color, orig_color = cg$color, self = as.numeric(vertex_vector),
    ideal = as.numeric(ideal_vector), hidden = !show,
    title = paste0("<p>", poles, "<br>Self: ", round(vertex_vector,2), "</p>"),
    stringsAsFactors = FALSE
  )

  edges_raw <- .extract_edges(wmatrix)
  max_w <- if(nrow(edges_raw)>0) max(abs(edges_raw$weight)) else 1
  
  if (nrow(edges_raw) > 0) {
    edges <- data.frame(
      from = as.character(edges_raw$from), to = as.character(edges_raw$to),
      width = 2 * abs(edges_raw$weight), arrows = "to",
      color = ifelse(edges_raw$weight > 0, "grey", "#CD5C5C"),
      weight = edges_raw$weight, hidden = FALSE, stringsAsFactors = FALSE
    )
  } else {
    edges <- data.frame(from=character(0), to=character(0), stringsAsFactors=FALSE)
  }

  .get_all_layouts <- function(wm, v, aa) {
    ig <- igraph::graph_from_adjacency_matrix(wm, weight=TRUE, mode="directed")
    .sc <- function(m) { if(nrow(m)==0) return(m); for(i in 1:2){ r=range(m[,i]); s=r[2]-r[1]; if(s<1e-9) s=1; m[,i]=(m[,i]-r[1])/s - 0.5 }; m*850 }
    nids <- v$id
    list(
      "graphopt" = data.frame(id=nids, x=.sc(igraph::layout_with_graphopt(ig))[,1], y=.sc(igraph::layout_with_graphopt(ig))[,2]),
      "circle"   = data.frame(id=nids, x=.sc(igraph::layout_in_circle(ig))[,1], y=.sc(igraph::layout_in_circle(ig))[,2]),
      "mds"      = data.frame(id=nids, x=.sc(igraph::layout_with_mds(ig))[,1], y=.sc(igraph::layout_with_mds(ig))[,2]),
      "grid"     = data.frame(id=nids, x=.sc(igraph::layout_on_grid(ig))[,1], y=.sc(igraph::layout_on_grid(ig))[,2]),
      "tree"     = data.frame(id=nids, x=.sc(igraph::layout_as_tree(ig, circular=TRUE))[,1], y=.sc(igraph::layout_as_tree(ig, circular=TRUE))[,2])
    )
  }

  g <- visNetwork::visNetwork(vertex, edges, height = height, width = width) %>%
    visNetwork::visIgraphLayout(layout = ifelse(layout=="graphopt","layout_with_graphopt","layout_in_circle"), randomSeed = 33) %>%
    visNetwork::visOptions(highlightNearest = list(enabled=TRUE, degree=0), selectedBy = "group")

  if (interactive_options) {
    js_panel <- "
    function(el, x) {
      el.style.height = '90vh';
      var network = this.network;

      var createPanel = function(id, title, pos, width) {
        var p = document.createElement('div'); p.id = id;
        Object.assign(p.style, {
          position: 'absolute', zIndex: '1000', backgroundColor: 'rgba(255,255,255,0.95)',
          padding: '10px', borderRadius: '8px', boxShadow: '0 2px 15px rgba(0,0,0,0.15)',
          border: '1px solid #ddd', fontFamily: 'Segoe UI, sans-serif', fontSize: '11px',
          width: width || '200px', maxHeight: '40px', overflowY: 'auto', transition: 'all 0.3s ease'
        }, pos);
        var head = document.createElement('div');
        head.style = 'display:flex; justify-content:space-between; align-items:center; cursor:pointer; border-bottom:1px solid #eee; padding-bottom:5px; margin-bottom:10px;';
        head.innerHTML = '<b>' + title + '</b><span class=\"tgl\">+</span>';
        p.appendChild(head);
        var body = document.createElement('div'); body.style.display = 'none'; p.appendChild(body);
        head.onclick = function() {
          var show = body.style.display === 'none'; body.style.display = show ? 'block' : 'none';
          head.querySelector('.tgl').innerText = show ? '−' : '+'; p.style.maxHeight = show ? '45%' : '40px';
        };
        el.appendChild(p); return body;
      };

      var visBody = createPanel('p_vis', 'Visualization Options', {top:'10px', right:'10px'});
      var simBody = null, timeline = null;
      if (x.sim_data) {
        simBody = createPanel('p_sim', 'Sim Context', {bottom:'80px', right:'10px', width:'240px'});
        timeline = document.createElement('div');
        Object.assign(timeline.style, {
          position:'absolute', bottom:'10px', left:'50%', transform:'translateX(-50%)',
          width:'50%', minWidth:'350px', backgroundColor:'rgba(255,255,255,0.95)',
          padding:'12px 25px', borderRadius:'35px', boxShadow:'0 4px 15px rgba(0,0,0,0.1)',
          border:'1px solid #ddd', zIndex:'1001', display:'flex', alignItems:'center', gap:'15px',
          fontFamily: 'Segoe UI, sans-serif'
        });
        el.appendChild(timeline);
      }

      visBody.innerHTML = '<label>Palette</label><select id=\"ps\" style=\"width:100%;\"></select>' +
        '<label>Layout</label><select id=\"ls\" style=\"width:100%;\"></select>' +
        '<div style=\"margin-top:10px;\"><b>Edge Filter</b> <span id=\"wt\">0.00</span><input type=\"range\" id=\"ws\" min=\"0\" max=\"'+x.max_weight.toFixed(2)+'\" step=\"0.01\" value=\"0\" style=\"width:100%;\"></div>' +
        '<div id=\"nl\" style=\"max-height:100px; overflow-y:auto; border:1px solid #eee; margin-top:5px;\"></div>';

      var ps = visBody.querySelector('#ps'), ls = visBody.querySelector('#ls');
      Object.keys(x.color_palette_js).forEach(k => { var o = document.createElement('option'); o.value=k; o.text=k; if(k===x.initial_palette) o.selected=true; ps.appendChild(o); });
      Object.keys(x.layouts).forEach(k => { var o = document.createElement('option'); o.value=k; o.text=k; if(k===x.initial_layout) o.selected=true; ls.appendChild(o); });

      if (simBody) {
        timeline.innerHTML = '<button id=\"ply\" style=\"background:#3498db; color:white; border:none; border-radius:16px; width:32px; height:32px; cursor:pointer; flex-shrink:0;\">▶</button>' +
          '<button id=\"pse\" style=\"background:#eee; border:1px solid #ccc; border-radius:16px; width:32px; height:32px; cursor:pointer; flex-shrink:0;\">II</button>' +
          '<div style=\"flex:1; display:flex; flex-direction:column;\"><div style=\"display:flex; justify-content:space-between; font-size:10px; color:#666;\"><b>Timeline</b> <b id=\"il\">Iter 0</b></div><input type=\"range\" id=\"sl\" min=\"0\" max=\"'+x.sim_data.max_iter+'\" value=\"0\" style=\"width:100%; accent-color:#3498db;\"></div>' +
          '<button id=\"rst\" style=\"background:white; border:1px solid #ddd; padding:4px 8px; border-radius:12px; font-size:10px;\">Reset</button>';
        simBody.innerHTML = '<b>Target Baseline</b><div id=\"al\" style=\"max-height:200px; overflow-y:auto;\"></div>';
        var sim = x.sim_data, hist = [], curI = 0, targetS = [...sim.initial_self];
        sim.lpoles.forEach((lp, i) => {
          var d = document.createElement('div'); d.style.marginBottom = '5px';
          d.innerHTML = '<div style=\"display:flex; justify-content:space-between; font-size:9px;\"><span>'+lp+'</span><b class=\"vb\">'+sim.initial_self[i].toFixed(2)+'</b></div><input type=\"range\" class=\"ts\" data-idx=\"'+i+'\" min=\"-1\" max=\"1\" step=\"0.05\" value=\"'+sim.initial_self[i]+'\" style=\"width:100%; height:4px;\">';
          simBody.querySelector('#al').appendChild(d);
        });
        var run = function() {
          var s = [...sim.initial_self], a = targetS.map((t, idx) => t - s[idx]), h = [[...s]], m = sim.max_iter;
          var thr = (v) => sim.threshold==='tanh'?Math.tanh(v):(sim.threshold==='saturation'?Math.max(-1,Math.min(1,v)):v);
          for(var i=0; i<m; i++){
            var next = s.map((v,idx) => thr(v+a[idx])), delta = next.map((v,idx) => v-s[idx]); s=next; var nA = new Array(s.length).fill(0);
            for(var r=0; r<s.length; r++) for(var c=0; c<s.length; c++) nA[r]+=sim.weights[c][r]*delta[c]; a=nA; h.push([...s]);
          }
          hist = h; update(curI);
        };
        var update = function(idx) {
          curI = idx; var data = hist[idx]; if(!data) return;
          timeline.querySelector('#sl').value = idx; timeline.querySelector('#il').innerText = 'Iter ' + idx;
          var pal = x.color_palette_js[ps.value];
          network.body.data.nodes.update(network.body.data.nodes.get().map((node, i) => {
            var v = data[i], idl = sim.initial_ideal[i], color = (idl===0)?pal[3]:(v===0?pal[2]:(Math.sign(v)===Math.sign(idl)?pal[1]:pal[0]));
            var pole = (v<0)?sim.lpoles[i]:(v>0?sim.rpoles[i]:sim.lpoles[i]+' - '+sim.rpoles[i]), sz = 20 + (40 * Math.abs(v));
            return { id:node.id, color:{background:color, border:'#222'}, label: pole+'\n('+v.toFixed(2)+')', size:sz, font: { vadjust: -sz - 15, strokeWidth: 4, strokeColor: '#fff' } };
          }));
        };
        timeline.querySelector('#sl').oninput = function() { update(parseInt(this.value)); };
        var timer = null;
        timeline.querySelector('#ply').onclick = function() { if(timer) clearInterval(timer); timer = setInterval(() => { if(curI < sim.max_iter) update(curI+1); else clearInterval(timer); }, 600); };
        timeline.querySelector('#pse').onclick = function() { clearInterval(timer); };
        timeline.querySelector('#rst').onclick = function() { targetS = [...sim.initial_self]; simBody.querySelectorAll('.ts').forEach((s, idx) => { s.value = targetS[idx]; s.parentNode.querySelector('.vb').innerText = targetS[idx].toFixed(2); }); run(); update(0); };
        simBody.querySelector('#al').oninput = function(e) { if(e.target.classList.contains('ts')) { var i = parseInt(e.target.dataset.idx); targetS[i] = parseFloat(e.target.value); e.target.parentNode.querySelector('.vb').innerText = targetS[i].toFixed(2); run(); } };
        run();
      }

      var upPal = function(m) {
        var pal = x.color_palette_js[m], nodesDS = network.body.data.nodes;
        nodesDS.update(nodesDS.get().map(n => {
          var idx = parseInt(n.id)-1, v = (simBody && hist[curI]) ? hist[curI][idx] : n.self, idl = (simBody) ? x.sim_data.initial_ideal[idx] : n.ideal;
          return { id:n.id, color:{background:(idl===0?pal[3]:(v===0?pal[2]:(Math.sign(v)===Math.sign(idl)?pal[1]:pal[0])))}};
        }));
      };
      ps.onchange = function() { upPal(this.value); };
      ls.onchange = function() {
        var l = x.layouts[this.value]; if(!l) return;
        network.setOptions({physics:{enabled:false}}); network.body.data.nodes.update(l.id.map((id,i)=>({id:String(id), x:l.x[i], y:l.y[i]})));
        upPal(ps.value); network.fit({animation:true});
      };
      visBody.querySelector('#ws').oninput = function() {
        var th = parseFloat(this.value); visBody.querySelector('#wt').innerText = th.toFixed(2);
        network.body.data.edges.update(network.body.data.edges.get().map(e => ({id:e.id, hidden:Math.abs(e.weight)<th})));
      };
      
      var nodesDS = network.body.data.nodes;
      nodesDS.get().forEach(n => {
        var d = document.createElement('div'); d.style.fontSize='10px';
        d.innerHTML = '<label><input type=\"checkbox\" class=\"nc\" data-id=\"'+n.id+'\" '+(n.hidden?'':'checked')+'> '+n.label.split('\n')[0]+'</label>';
        visBody.querySelector('#nl').appendChild(d);
      });
      visBody.querySelector('#nl').onchange = function(e) { if(e.target.classList.contains('nc')) { var n = nodesDS.get(e.target.dataset.id); n.hidden=!e.target.checked; nodesDS.update(n); } };
      
      var btn = document.createElement('button'); btn.innerHTML = '📷 Export PNG';
      Object.assign(btn.style, { position:'absolute', bottom:'15px', left:'15px', zIndex:'1000', padding:'8px 15px', background:'#fff', border:'1px solid #ddd', borderRadius:'20px', cursor:'pointer', fontSize:'11px' });
      btn.onclick = function() { var c = el.getElementsByTagName('canvas')[0], l = document.createElement('a'); l.download = 'WimpTools.png'; l.href = c.toDataURL('image/png', 1.0); l.click(); };
      el.appendChild(btn);
    }
    "
    g$x$interactive_options <- interactive_options
    g$x$min_weight    <- min_weight
    g$x$max_weight    <- max_w
    g$x$layouts       <- .get_all_layouts(wd$wmatrix, vertex, area_attr)
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
#' @description Plot the ideal self based on the constructs and their relations.
#' @param wimp Subject's WimpGrid
#' @param ... Additional arguments
#' @export
idealdigraph <- function(wimp, ...) {
  digraph(wimp = wimp, vertex_vector = wimp$vertices$ideal, ...)
}

# Simulation Digraph ---------------------------------------------------------

#' Simulation digraph -- simdigraph()
#' @description Plot the hypothetical self for a given scenario iteration.
#' @param scn A scenario matrix or wimp object.
#' @param ... Additional arguments
#' @export
simdigraph <- function(scn, ...) {
  if (inherits(scn, "wimp")) {
    wimp <- .align_wimp(scn, exclude_dilemmatics = FALSE)
    sim_data <- list(initial_self = wimp$vertices$self, initial_ideal = wimp$vertices$ideal, weights = as.matrix(wimp$global$weight_matrix), threshold = "saturation", max_iter = 10, lpoles = wimp$vertices$left_pole, rpoles = wimp$vertices$right_pole)
    return(digraph(wimp = wimp, vertex_vector = wimp$vertices$self, sim_data = sim_data, ...))
  } else if (inherits(scn, "scn")) {
    sim_data <- list(initial_self = scn$self$self, initial_ideal = scn$self$ideal, weights = as.matrix(scn$weights), threshold = scn$params$threshold, max_iter = scn$params$max_iter, lpoles = scn$constructs$left_pole, rpoles = scn$constructs$right_pole)
    return(digraph(wimp = scn, vertex_vector = scn$values[1, ], sim_data = sim_data, ...))
  }
}
