# Shared responsive CSS for widget chrome (tabs, floating icon buttons, tables).
# Every widget renders inside its own iframe, so `vw`/`vh` units and `@media`
# queries below are scoped to that widget's own rendered size, not the browser
# window - shrinking the GridStack card shrinks these exactly like a container
# query would, without needing @container support.
.wt_responsive_css <- function() "
  .wt-tab-header button, .wct-tab-header button, .wsim-tab-btn, .wrg-tab-header button, .wt-tab-btn {
    white-space: nowrap;
    overflow: hidden;
    text-overflow: ellipsis;
    min-width: 0;
  }
  table { font-size: 13px; }
  @media (max-width: 520px) {
    .wt-tab-header button, .wct-tab-header button, .wsim-tab-btn, .wrg-tab-header button, .wt-tab-btn {
      padding: 8px 6px !important;
      font-size: 12px !important;
    }
    table { font-size: 12px; }
    table th, table td { padding: 7px 8px !important; }
    /* Leave a gutter for the floating info/export/fullscreen buttons so they
       don't sit on top of plot traces or axis tick labels drawn near the
       right edge (e.g. the RepGrid dilemmas dumbbell chart). Deliberately
       excludes .wt-tab-content > .html-widget: in widget_digraph the
       floating buttons are appended INSIDE that same html-widget div, so
       shrinking it drags the buttons in with it, leaving a dead gap between
       them and the widget's actual right edge instead of a gutter. */
    .wct-tab-content > .html-widget,
    .wsim-tab-content .html-widget, .wrg-tab-content > .html-widget {
      width: calc(100% - 58px) !important;
    }
  }
  @media (max-width: 340px) {
    .wt-tab-header button, .wct-tab-header button, .wsim-tab-btn, .wrg-tab-header button, .wt-tab-btn {
      padding: 6px 4px !important;
      font-size: 10.5px !important;
    }
    table { font-size: 11px; }
    table th, table td { padding: 5px 6px !important; }
    .wct-tab-content > .html-widget,
    .wsim-tab-content .html-widget, .wrg-tab-content > .html-widget {
      width: calc(100% - 48px) !important;
    }
  }
"

#' Generate Tabbed HTML Widget for Psychlab
#'
#' @description
#' Creates a standalone HTML widget containing a two-tab interface for Psychlab integration.
#' The first tab displays the interactive digraph network, and the second tab displays
#' the weight matrix heatmap.
#'
#' @param x A \code{wimp} object created by \code{\link{wimp}}.
#' @param ... Additional arguments passed to \code{digraph}.
#'
#' @return A `browsable` HTML object containing the tabbed interface.
#' @export
#'
#' @importFrom htmltools tagList tags browsable HTML
#'
#' @examples
#' \dontrun{
#' # Assuming 'w' is a valid wimp object
#' widget <- widget_digraph(w)
#' widget
#' }
widget_digraph <- function(x, lang = "en", ...) {
  if (!inherits(x, "wimp")) stop("Input must be a 'wimp' object.")
  if (!lang %in% c("en", "es")) lang <- "en"
  
  export_name <- deparse(substitute(x))
  
  # --- Localization (see R/i18n.R) ---
  t <- wt_i18n(lang)

  # 1. Generate individual widgets
  g <- digraph(x, export_name = export_name, lang = lang, ...)
  h <- weight_heatmap(x, lang = lang)
  
  # 2. Define custom CSS for tabs
  css <- "
    body, html {
      margin: 0;
      padding: 0;
      width: 100%;
      height: 100%;
      overflow: hidden;
    }
    * {
      box-sizing: border-box;
    }
    :root {
      accent-color: #8cc63f !important;
    }
    input[type=\"checkbox\"], input[type=\"radio\"], input[type=\"range\"] {
      accent-color: #8cc63f !important;
    }
    select:focus { 
      border-color: #8cc63f !important; 
      outline: none; 
    }
    option:checked { 
      background-color: #8cc63f !important; 
      color: white !important; 
    }
    option:hover, option:focus, option:active { 
      background-color: #8cc63f !important; 
      color: white !important;
      box-shadow: 0 0 10px 100px #8cc63f inset !important;
    }
    ::selection { 
      background-color: #8cc63f !important; 
      color: white !important; 
    }
    .wt-tab-container {
      width: 100%;
      height: 100vh;
      display: flex;
      flex-direction: column;
      font-family: 'Inter', 'Roboto', 'Segoe UI', sans-serif;
    }
    .wt-tab-header {
      overflow: hidden;
      border: 1px solid #eaeaea;
      background-color: #ffffff;
      display: flex;
      border-radius: 4px 4px 0 0;
    }
    .wt-tab-header button {
      background-color: inherit;
      border: none;
      outline: none;
      cursor: pointer;
      padding: 12px 24px;
      transition: 0.3s;
      font-size: 14px;
      font-weight: 600;
      color: #666666;
      flex-grow: 1;
    }
    .wt-tab-header button:hover {
      background-color: #f9f9f9;
    }
    .wt-tab-header button.active {
      background-color: #ffffff;
      color: #8cc63f;
      border-bottom: 2px solid #8cc63f;
    }
    .wt-tab-content {
      position: relative;
      display: none;
      padding: 0;
      border: 1px solid #ccc;
      border-top: none;
      flex-grow: 1;
      height: 0;
      background-color: #ffffff;
      border-radius: 0 0 4px 4px;
    }
    .wt-tab-content > .html-widget {
      width: 100% !important;
      height: 100% !important;
      flex-grow: 1;
    }
  "
  css <- paste0(css, .wt_responsive_css())

  # 3. Define vanilla JS for tab switching
  js <- "
    function openPsychlabTab(evt, tabName) {
      var i, tabcontent, tablinks;
      
      // Hide all tab content
      tabcontent = document.getElementsByClassName('wt-tab-content');
      for (i = 0; i < tabcontent.length; i++) {
        tabcontent[i].style.display = 'none';
      }
      
      // Remove 'active' class from all buttons
      tablinks = document.getElementsByClassName('wt-tab-btn');
      for (i = 0; i < tablinks.length; i++) {
        tablinks[i].className = tablinks[i].className.replace(' active', '');
      }
      
      // Show the current tab, and add an 'active' class to the button that opened the tab
      document.getElementById(tabName).style.display = 'flex';
      document.getElementById(tabName).style.flexDirection = 'column';
      evt.currentTarget.className += ' active';
      
      // Trigger resize event so htmlwidgets (plotly, visNetwork) render correctly in newly visible containers
      window.dispatchEvent(new Event('resize'));
    }
    
    function downloadHeatmap() {
      var plotEl = document.querySelector('#psychlab_tab_heatmap .js-plotly-plot');
      if (plotEl && window.Plotly) {
        Plotly.downloadImage(plotEl, {format: 'png', width: plotEl.clientWidth, height: plotEl.clientHeight, filename: 'WIMP_EXPORT_NAME_Weight_Matrix'});
      }
    }
  "
  js <- gsub("WIMP_EXPORT_NAME", export_name, js)
  
  # 4. Construct UI using htmltools
  ui <- htmltools::tagList(
    htmltools::tags$style(htmltools::HTML(css)),
    htmltools::tags$script(htmltools::HTML(js)),
    
    htmltools::tags$div(class = "wt-tab-container",
      htmltools::tags$div(class = "wt-tab-header",
        htmltools::tags$button(
          class = "wt-tab-btn active", 
          onclick = "openPsychlabTab(event, 'psychlab_tab_graph')", 
          t$self_digraph
        ),
        htmltools::tags$button(
          class = "wt-tab-btn", 
          onclick = "openPsychlabTab(event, 'psychlab_tab_heatmap')", 
          t$weight_matrix
        )
      ),
      
      htmltools::tags$div(id = "psychlab_tab_graph", class = "wt-tab-content", style = "display:flex; flex-direction:column;", 
        g
      ),
      htmltools::tags$div(id = "psychlab_tab_heatmap", class = "wt-tab-content", 
        h,
        htmltools::HTML(paste0("
          <div id='hm_btn_container' style='position:absolute;bottom:15px;right:15px;display:flex;flex-direction:column;gap:8px;z-index:1000;align-items:center;'>
            <div style='background-color: rgba(255, 255, 255, 0.95); width: clamp(24px, 4vmin, 32px); height: clamp(24px, 4vmin, 32px); border-radius: 6px; box-shadow: 0 2px 10px rgba(0,0,0,0.1); border: 1px solid #ddd; display: flex; align-items: center; justify-content: center; cursor: pointer; transition: all 0.2s;' title='", t$info, "' onmouseover=\"this.style.backgroundColor='#f5f5f5'\" onmouseout=\"this.style.backgroundColor='rgba(255, 255, 255, 0.95)'\" onclick=\"var m=document.getElementById('heatmap_info_modal'); m.style.display=(m.style.display==='block'?'none':'block');\">
              <svg width='60%' height='60%' viewBox='0 0 24 24' fill='none' stroke='#333' stroke-width='2.5' stroke-linecap='round' stroke-linejoin='round'><circle cx='12' cy='12' r='10'></circle><line x1='12' y1='16' x2='12' y2='12'></line><line x1='12' y1='8' x2='12.01' y2='8'></line></svg>
            </div>
            <div id='hm_settings_placeholder'></div>
            <div style='background-color: rgba(255, 255, 255, 0.95); width: clamp(24px, 4vmin, 32px); height: clamp(24px, 4vmin, 32px); border-radius: 6px; box-shadow: 0 2px 10px rgba(0,0,0,0.1); border: 1px solid #ddd; display: flex; align-items: center; justify-content: center; cursor: pointer; transition: all 0.2s;' title='", t$export_png, "' onmouseover=\"this.style.backgroundColor='#f5f5f5'\" onmouseout=\"this.style.backgroundColor='rgba(255, 255, 255, 0.95)'\" onclick=\"downloadHeatmap()\">
              <svg width='60%' height='60%' viewBox='0 0 24 24' fill='none' stroke='#333' stroke-width='2.5' stroke-linecap='round' stroke-linejoin='round'><path d='M21 15v4a2 2 0 0 1-2 2H5a2 2 0 0 1-2-2v-4'></path><polyline points='7 10 12 15 17 10'></polyline><line x1='12' y1='15' x2='12' y2='3'></line></svg>
            </div>
            <div style='background-color: rgba(255, 255, 255, 0.95); width: clamp(24px, 4vmin, 32px); height: clamp(24px, 4vmin, 32px); border-radius: 6px; box-shadow: 0 2px 10px rgba(0,0,0,0.1); border: 1px solid #ddd; display: flex; align-items: center; justify-content: center; cursor: pointer; transition: all 0.2s;' title='", t$fullscreen, "' onmouseover=\"this.style.backgroundColor='#f5f5f5'\" onmouseout=\"this.style.backgroundColor='rgba(255, 255, 255, 0.95)'\" onclick=\"var el=this.closest('.wt-tab-container')||this.closest('.wt-tab-content'); if(!document.fullscreenElement){el.requestFullscreen().catch(e=>console.log(e))}else{document.exitFullscreen()}\">
              <svg width='60%' height='60%' viewBox='0 0 24 24' fill='none' stroke='#333' stroke-width='2.5' stroke-linecap='round' stroke-linejoin='round'><path d='M8 3H5a2 2 0 0 0-2 2v3m18 0V5a2 2 0 0 0-2-2h-3m0 18h3a2 2 0 0 0 2-2v-3M3 16v3a2 2 0 0 0 2 2h3'></path></svg>
            </div>
          </div>
          
          <div id='heatmap_info_modal' style='position: absolute; top: 50%; left: 50%; transform: translate(-50%, -50%); width: 90%; max-width: 450px; max-height: 80vh; overflow-y: auto; box-sizing: border-box; background-color: #fff; z-index: 2000; padding: 20px; border-radius: 8px; box-shadow: 0 4px 20px rgba(0,0,0,0.2); border: 1px solid #eaeaea; display: none; font-family: Inter, Roboto, sans-serif;'>
            <div style='display:flex; justify-content:space-between; align-items:center; border-bottom:1px solid #eaeaea; padding-bottom:10px; margin-bottom:15px;'>
              <h3 style='margin:0; color:#444; font-size:16px;'>", t$weight_matrix, "</h3>
              <span style='cursor:pointer; font-size:20px; font-weight:bold; color:#888; line-height:1;' onclick=\"document.getElementById('heatmap_info_modal').style.display='none';\">&times;</span>
            </div>
            <p style='margin:0; color:#666; font-size:13px; line-height:1.6;'>", t$info_text_heatmap, "</p>
          </div>
        "))
      )
    )
  )
  
  # 5. Make it browsable so it renders in RStudio viewer or standalone HTML
  return(htmltools::browsable(ui))
}


# widget_centrality -----------------------------------------------------------

#' Centrality Widget for Psychlab
#'
#' @description Creates a two-tab HTML widget: Tab 1 shows the PB Plot
#'   (Presence-Balance space), Tab 2 shows a sortable table with P, B, Degree,
#'   Closeness and Betweenness for each construct.
#'
#' @param x A \code{wimp} object.
#' @param lang Language for the UI. \code{"en"} (default) or \code{"es"}.
#' @param ... Additional arguments passed to \code{pb_plot}.
#'
#' @return A \code{browsable} HTML object with a two-tab interface.
#' @export
#'
#' @importFrom htmltools tagList tags browsable HTML
#'
#' @examples
#' \dontrun{
#' widget_centrality(example_wimp, lang = "es")
#' }
widget_centrality <- function(x, lang = "en", ...) {
  if (!inherits(x, "wimp")) stop("Input must be a 'wimp' object.")
  if (!lang %in% c("en", "es")) lang <- "en"

  t <- wt_i18n(lang)

  # 1. Compute metrics --------------------------------------------------------
  pb  <- as.data.frame(pb_index(x))
  deg <- as.data.frame(degree_index(x, method = "wnorm"))
  cl  <- as.data.frame(close_index(x))
  bw  <- as.data.frame(betw_index(x))

  tbl <- data.frame(
    construct   = rownames(pb),
    p           = round(pb$p,           3),
    b           = round(pb$b,           3),
    degree      = round(deg$All,        3),
    closeness   = round(cl$Closeness,   3),
    betweenness = round(bw$Betweenness, 3),
    stringsAsFactors = FALSE
  )

  # 2. Build HTML table -------------------------------------------------------
  col_headers <- c(t$construct, t$presence, t$balance,
                   t$degree, t$closeness, t$betweenness)
  col_keys    <- c("construct", "p", "b", "degree", "closeness", "betweenness")

  header_cells <- paste(
    mapply(function(label, key) {
      sprintf(
        "<th onclick=\"sortWctTable('%s')\" style='cursor:pointer;user-select:none;'>%s <span id='wct_sort_%s' style='font-size:15px;color:#aaa;'>&#8597;</span></th>",
        key, label, key)
    }, col_headers, col_keys),
    collapse = "\n")

  data_rows <- paste(apply(tbl, 1, function(row) {
    sprintf(
      "<tr><td>%s</td><td style='text-align:right;'>%s</td><td style='text-align:right;'>%s</td><td style='text-align:right;'>%s</td><td style='text-align:right;'>%s</td><td style='text-align:right;'>%s</td></tr>",
      row["construct"], row["p"], row["b"],
      row["degree"], row["closeness"], row["betweenness"])
  }), collapse = "\n")

  table_html <- paste0(
    "<div style='width:100%; height:100%; display:flex; flex-direction:column; font-family:Inter,Roboto,sans-serif;'>",
    "<div style='padding:16px 20px 8px 20px; flex-shrink:0;'>",
    "<input id='wct_search' type='text' placeholder='", t$construct, "...'",
    " oninput='filterWctTable()'",
    " style='width:100%; max-width:400px; padding:8px 12px; border:1px solid #ddd;",
    " border-radius:6px; font-size:13px; outline:none; font-family:Inter,Roboto,sans-serif;'>",
    "</div>",
    "<div style='flex:1; overflow:auto; padding:0 clamp(16px, 10vw, 70px) 60px 20px;'>",
    "<table id='wct_table' style='width:100%; border-collapse:collapse; font-size:13px;'>",
    "<thead><tr style='background:#f8f9fa; border-bottom:2px solid #8cc63f;'>",
    header_cells,
    "</tr></thead>",
    "<tbody id='wct_tbody'>", data_rows, "</tbody>",
    "</table></div></div>"
  )

  # 3. PB Plot ----------------------------------------------------------------
  pb_p <- pb_plot(x, lang = lang, ...) %>%
    plotly::config(displayModeBar = FALSE)

  # 4. CSS (unique prefix wct- to avoid conflicts with widget_digraph) --------
  css <- "
    body, html { margin: 0; padding: 0; width: 100%; height: 100%; overflow: hidden; }
    * { box-sizing: border-box; }
    .wct-container {
      width: 100%;
      height: 100%;
      display: flex;
      flex-direction: column;
      font-family: 'Inter', 'Roboto', 'Segoe UI', sans-serif;
    }
    .wct-tab-header {
      overflow: hidden;
      border: 1px solid #eaeaea;
      background-color: #ffffff;
      display: flex;
      border-radius: 4px 4px 0 0;
    }
    .wct-tab-header button {
      background-color: inherit;
      border: none;
      outline: none;
      cursor: pointer;
      padding: 12px 24px;
      transition: 0.3s;
      font-size: 14px;
      font-weight: 600;
      color: #666666;
      flex-grow: 1;
    }
    .wct-tab-header button:hover { background-color: #f9f9f9; }
    .wct-tab-header button.wct-active {
      background-color: #ffffff;
      color: #8cc63f;
      border-bottom: 2px solid #8cc63f;
    }
    .wct-tab-content {
      position: relative;
      display: none;
      padding: 0;
      border: 1px solid #ccc;
      border-top: none;
      flex-grow: 1;
      height: 0;
      background-color: #ffffff;
      border-radius: 0 0 4px 4px;
    }
    .wct-tab-content > .html-widget {
      width: 100% !important;
      height: 100% !important;
      flex-grow: 1;
    }
    #wct_table thead th { padding:10px 14px; text-align:left; font-weight:600; color:#555; font-size:12px; white-space:nowrap; }
    #wct_table thead th:hover { background:#eef7e0; }
    #wct_table tbody tr { border-bottom:1px solid #f0f0f0; }
    #wct_table tbody tr:hover { background:#f4f9ef; }
    #wct_table tbody td { padding:9px 14px; color:#333; }
    #wct_search:focus { border-color:#8cc63f !important; box-shadow:0 0 0 2px rgba(140,198,63,0.2); }
  "
  css <- paste0(css, .wt_responsive_css())

  # 5. JS --------------------------------------------------------------------
  js <- sprintf("
    function openWctTab(evt, tabName) {
      var i, tabcontent, tablinks;
      tabcontent = document.getElementsByClassName('wct-tab-content');
      for (i = 0; i < tabcontent.length; i++) { tabcontent[i].style.display = 'none'; }
      tablinks = document.getElementsByClassName('wct-tab-btn');
      for (i = 0; i < tablinks.length; i++) {
        tablinks[i].className = tablinks[i].className.replace(' wct-active', '');
      }
      document.getElementById(tabName).style.display = 'flex';
      document.getElementById(tabName).style.flexDirection = 'column';
      evt.currentTarget.className += ' wct-active';
      window.dispatchEvent(new Event('resize'));
    }

    function downloadPBPlot() {
      var plotEl = document.querySelector('#wct_tab_pb .js-plotly-plot');
      if (plotEl && window.Plotly) {
        Plotly.downloadImage(plotEl, {format: 'png', width: plotEl.clientWidth, height: plotEl.clientHeight, filename: 'WimpTools_PB_Plot'});
      }
    }

    function downloadCentralityCSV() {
      var rows = [['%s','%s','%s','%s','%s','%s']];
      document.querySelectorAll('#wct_tbody tr').forEach(function(r) {
        if (r.style.display !== 'none')
          rows.push(Array.from(r.cells).map(function(c){ return c.innerText.trim(); }));
      });
      var csv = rows.map(function(r){ return r.join(','); }).join('\\n');
      var a = document.createElement('a');
      a.href = 'data:text/csv;charset=utf-8,' + encodeURIComponent(csv);
      a.download = 'WimpTools_Centrality.csv';
      a.click();
    }

    var _wctDir = {};
    function sortWctTable(col) {
      var tbody = document.getElementById('wct_tbody');
      var rows  = Array.from(tbody.querySelectorAll('tr'));
      var cols  = ['construct','p','b','degree','closeness','betweenness'];
      var idx   = cols.indexOf(col); if (idx < 0) return;
      _wctDir[col] = !_wctDir[col]; var asc = _wctDir[col];
      rows.sort(function(a, b) {
        var va = a.cells[idx].innerText.trim(), vb = b.cells[idx].innerText.trim();
        var na = parseFloat(va), nb = parseFloat(vb);
        if (!isNaN(na) && !isNaN(nb)) return asc ? na - nb : nb - na;
        return asc ? va.localeCompare(vb) : vb.localeCompare(va);
      });
      document.querySelectorAll('[id^=wct_sort_]').forEach(function(e){ e.innerHTML='&#8597;'; e.style.color='#aaa'; });
      document.getElementById('wct_sort_'+col).innerHTML = asc ? '&#8593;' : '&#8595;';
      document.getElementById('wct_sort_'+col).style.color = '#8cc63f';
      rows.forEach(function(r){ tbody.appendChild(r); });
    }

    function filterWctTable() {
      var q = document.getElementById('wct_search').value.toLowerCase();
      document.querySelectorAll('#wct_tbody tr').forEach(function(r) {
        r.style.display = r.innerText.toLowerCase().includes(q) ? '' : 'none';
      });
    }
  ",
  t$construct, t$presence, t$balance, t$degree, t$closeness, t$betweenness)

  # 6. Helper builders -------------------------------------------------------
  svg_dl <- "<path d='M21 15v4a2 2 0 0 1-2 2H5a2 2 0 0 1-2-2v-4'></path><polyline points='7 10 12 15 17 10'></polyline><line x1='12' y1='15' x2='12' y2='3'></line>"
  svg_i  <- "<circle cx='12' cy='12' r='10'></circle><line x1='12' y1='16' x2='12' y2='12'></line><line x1='12' y1='8' x2='12.01' y2='8'></line>"
  svg_fs <- "<path d='M8 3H5a2 2 0 0 0-2 2v3m18 0V5a2 2 0 0 0-2-2h-3m0 18h3a2 2 0 0 0 2-2v-3M3 16v3a2 2 0 0 0 2 2h3'></path>"

  .btn <- function(right_px, title_txt, onclick_fn, svg_body) {
    htmltools::HTML(paste0(
      "<div style='background-color:rgba(255,255,255,0.95);width:clamp(26px, 4vmin, 34px);height:clamp(26px, 4vmin, 34px);border-radius:6px;",
      "box-shadow:0 2px 10px rgba(0,0,0,0.1);border:1px solid #ddd;display:flex;",
      "align-items:center;justify-content:center;cursor:pointer;transition:all 0.2s;'",
      " title='", title_txt, "'",
      " onmouseover=\"this.style.backgroundColor='#f5f5f5'\"",
      " onmouseout=\"this.style.backgroundColor='rgba(255,255,255,0.95)'\"",
      " onclick=\"", onclick_fn, "\">",
      "<svg width='60%' height='60%' viewBox='0 0 24 24' fill='none' stroke='#333'",
      " stroke-width='2.5' stroke-linecap='round' stroke-linejoin='round'>",
      svg_body, "</svg></div>"
    ))
  }

  .modal <- function(id, title_txt, body_txt) {
    htmltools::HTML(paste0(
      "<div id='", id, "' style='position:absolute;top:50%;left:50%;",
      "transform:translate(-50%,-50%);width:90%;max-width:450px;max-height:80vh;overflow-y:auto;box-sizing:border-box;background:#fff;",
      "z-index:2000;padding:20px;border-radius:8px;",
      "box-shadow:0 4px 20px rgba(0,0,0,0.2);border:1px solid #eaeaea;display:none;",
      "font-family:Inter,Roboto,sans-serif;'>",
      "<div style='display:flex;justify-content:space-between;align-items:center;",
      "border-bottom:1px solid #eaeaea;padding-bottom:10px;margin-bottom:15px;'>",
      "<h3 style='margin:0;color:#444;font-size:16px;'>", title_txt, "</h3>",
      "<span style='cursor:pointer;font-size:20px;font-weight:bold;color:#888;line-height:1;'",
      " onclick=\"document.getElementById('", id, "').style.display='none';\">&times;</span>",
      "</div>",
      "<p style='margin:0;color:#666;font-size:13px;line-height:1.6;'>", body_txt, "</p>",
      "</div>"
    ))
  }

  fs_onclick <- "var el=this.closest('.wct-container'); if(!document.fullscreenElement){el.requestFullscreen().catch(e=>console.log(e))}else{document.exitFullscreen()}"

  # 7. Assemble UI -----------------------------------------------------------
  ui <- htmltools::tagList(
    htmltools::tags$style(htmltools::HTML(css)),
    htmltools::tags$script(htmltools::HTML(js)),

    htmltools::tags$div(class = "wct-container",

      htmltools::tags$div(class = "wct-tab-header",
        htmltools::tags$button(
          class   = "wct-tab-btn wct-active",
          onclick = "openWctTab(event, 'wct_tab_pb')",
          t$pb_plot_tab
        ),
        htmltools::tags$button(
          class   = "wct-tab-btn",
          onclick = "openWctTab(event, 'wct_tab_table')",
          t$centrality_tab
        )
      ),

      # Tab 1 - PB Plot
      htmltools::tags$div(id = "wct_tab_pb", class = "wct-tab-content",
        style = "display:flex; flex-direction:column;",
        pb_p,
        .modal("pb_info_modal", t$pb_plot_tab, t$info_text_pb),
        htmltools::HTML("<div style='position:absolute;bottom:15px;right:15px;z-index:1000;display:flex;flex-direction:column;gap:8px;align-items:center;'>"),
        .btn(0, t$info,       "var m=document.getElementById('pb_info_modal'); m.style.display=(m.style.display==='block'?'none':'block');", svg_i),
        .btn(0, t$export_png, "downloadPBPlot()",                                                svg_dl),
        .btn(0, t$fullscreen, fs_onclick,                                                        svg_fs),
        htmltools::HTML("</div>")
      ),

      # Tab 2 - Centrality table
      htmltools::tags$div(id = "wct_tab_table", class = "wct-tab-content",
        htmltools::HTML(table_html),
        .modal("ct_info_modal", t$centrality_tab, t$info_text_centrality),
        htmltools::HTML("<div style='position:absolute;bottom:15px;right:15px;z-index:1000;display:flex;flex-direction:column;gap:8px;align-items:center;'>"),
        .btn(0, t$info,       "var m=document.getElementById('ct_info_modal'); m.style.display=(m.style.display==='block'?'none':'block');", svg_i),
        .btn(0, t$export_png, "downloadCentralityCSV()",                                         svg_dl),
        .btn(0, t$fullscreen, fs_onclick,                                                        svg_fs),
        htmltools::HTML("</div>")
      )
    )
  )

  return(htmltools::browsable(ui))
}

# widget_simulation -----------------------------------------------------------

#' Simulation Widget for Psychlab
#'
#' @description Creates a two-tab HTML widget with a fixed sidebar for simulation controls.
#'   Tab 1 shows the Network View (simdigraph).
#'   Tab 2 shows the PCSD (Personal Construct System Dynamics) Plotly Chart.
#'
#' @param scn A \code{scn} or \code{wimp} object containing the scenario data.
#' @param lang Language for the UI. \code{"en"} (default) or \code{"es"}.
#' @param ... Additional arguments passed to \code{simdigraph}.
#'
#' @return A \code{browsable} HTML object with a two-tab interface and a fixed sidebar.
#' @export
#'
#' @importFrom htmltools tagList tags browsable HTML
#' @importFrom plotly plot_ly config layout
#'
#' @examples
#' \dontrun{
#' widget_simulation(example_wimp, lang = "es")
#' }
widget_simulation <- function(scn, lang = "en", ...) {
  if (!inherits(scn, c("wimp", "scn"))) stop("Input must be a 'wimp' or 'scn' object.")
  if (!lang %in% c("en", "es")) lang <- "en"

  t <- wt_i18n(lang)

  # 1. Main network plot (this will inject controls into #wsim_sidebar)
  network_plot <- simdigraph(scn, lang = lang, ...)

  # 2. Empty Plotly container for PCSD
  # We use plotly::plot_ly to ensure htmlwidgets loads the plotly library
  plotly_plot <- plotly::plot_ly(x = c(0), y = c(0), type = 'scatter', mode = 'lines', line = list(color = 'transparent'), showlegend = FALSE) %>%
    plotly::layout(
      title = list(text = "<b>PERSONAL CONSTRUCTS</b>"),
      xaxis = list(title = "ITERATIONS", zeroline = FALSE),
      yaxis = list(title = "SELF DIFFERENTIAL", zeroline = TRUE, zerolinecolor = "#666", zerolinewidth = 2, range = c(-2, 2))
    ) %>%
    plotly::config(displayModeBar = FALSE)

  # 3. CSS for layout
  css <- "
    body, html { margin: 0; padding: 0; width: 100%; height: 100%; overflow: hidden; }
    * { box-sizing: border-box; }
    .wsim-container {
      width: 100%;
      height: 100%;
      display: flex;
      flex-direction: row;
      font-family: 'Inter', 'Roboto', 'Segoe UI', sans-serif;
    }
    #wsim_sidebar {
      width: clamp(180px, 30vw, 320px);
      min-width: 160px;
      flex-shrink: 0;
      height: 100%;
      background: #f8f9fa;
      border-right: 1px solid #ddd;
      padding: clamp(8px, 2vw, 15px);
      overflow-y: auto;
      overflow-x: hidden;
      display: flex;
      flex-direction: column;
    }
    .wsim-main {
      flex: 1;
      min-width: 0;
      display: flex;
      flex-direction: column;
      height: 100%;
      overflow: hidden;
    }
    @media (max-width: 480px) {
      .wsim-container { flex-direction: column; }
      #wsim_sidebar {
        width: 100%;
        min-width: 0;
        height: auto;
        max-height: 45%;
        border-right: none;
        border-bottom: 1px solid #ddd;
      }
      .wsim-main { flex: 1; min-height: 0; }
    }
    .wsim-tab-header {
      overflow: hidden;
      border: 1px solid #eaeaea;
      background-color: #ffffff;
      display: flex;
      border-radius: 4px 4px 0 0;
    }
    .wsim-tab-btn {
      background-color: inherit;
      border: none;
      outline: none;
      cursor: pointer;
      padding: 12px 24px;
      transition: 0.3s;
      font-size: 14px;
      font-weight: 600;
      color: #666666;
      flex-grow: 1;
    }
    .wsim-tab-btn:hover {
      background-color: #f9f9f9;
    }
    .wsim-tab-btn.wsim-active {
      background-color: #ffffff;
      color: #8cc63f;
      border-bottom: 2px solid #8cc63f;
    }
    .wsim-tab-content {
      display: none;
      flex: 1;
      position: relative;
      min-height: 0;
    }
    .wsim-tab-content.wsim-active-tab {
      display: flex;
      flex-direction: column;
    }
    .wsim-tab-content .html-widget {
      flex: 1;
      min-height: 0;
      width: 100% !important;
      height: 100% !important;
    }
    .wsim-tab-content .html-widget .js-plotly-plot { height: 100% !important; }
  "
  css <- paste0(css, .wt_responsive_css())

  # 4. JS for tab switching
  js <- "
    function openWsimTab(evt, tabName) {
      document.querySelectorAll('.wsim-tab-content').forEach(function(e) {
        e.style.display = 'none';
        e.classList.remove('wsim-active-tab');
      });
      document.querySelectorAll('.wsim-tab-btn').forEach(function(b) {
        b.classList.remove('wsim-active');
      });
      var el = document.getElementById(tabName);
      el.style.display = 'flex';
      el.classList.add('wsim-active-tab');
      evt.currentTarget.classList.add('wsim-active');
      
      // Force plot resize if needed
      if (tabName === 'wsim_tab_pcsd') {
        setTimeout(function() {
          var pDiv = document.getElementById('wsim_pcsd_plot');
          if (pDiv && window.Plotly) Plotly.relayout(pDiv, {autosize: true});
        }, 50);
      } else if (tabName === 'wsim_tab_network') {
        window.dispatchEvent(new Event('resize'));
      }
    }
    
    function downloadSimNetwork() {
       // Fire a resize event first in case it's distorted
       window.dispatchEvent(new Event('resize'));
       setTimeout(function() {
          // Look for the visNetwork canvas
          var netCanvas = document.querySelector('#wsim_tab_network canvas');
          if (netCanvas) {
            var link = document.createElement('a');
            link.download = 'WIMP_EXPORT_NAME_Simulated_Network.png';
            link.href = netCanvas.toDataURL('image/png');
            link.click();
          }
       }, 500);
    }
    
    function downloadSimPcsd() {
      var plotEl = document.getElementById('wsim_pcsd_plot');
      if (plotEl && window.Plotly) {
        Plotly.downloadImage(plotEl, {format: 'png', width: plotEl.clientWidth, height: plotEl.clientHeight, filename: 'WIMP_EXPORT_NAME_Simulated_PCSD'});
      }
    }
  "

  svg_fs <- "<path d='M8 3H5a2 2 0 0 0-2 2v3m18 0V5a2 2 0 0 0-2-2h-3m0 18h3a2 2 0 0 0 2-2v-3M3 16v3a2 2 0 0 0 2 2h3'></path>"
  fs_onclick <- "var el=this.closest('.wsim-container'); if(!document.fullscreenElement){el.requestFullscreen().catch(e=>console.log(e))}else{document.exitFullscreen()}"

  .btn <- function(right_px, title_txt, onclick_fn, svg_body) {
    htmltools::HTML(paste0(
      "<div style='background-color:rgba(255,255,255,0.95);width:clamp(26px, 4vmin, 34px);height:clamp(26px, 4vmin, 34px);border-radius:6px;",
      "box-shadow:0 2px 10px rgba(0,0,0,0.1);border:1px solid #ddd;display:flex;",
      "align-items:center;justify-content:center;cursor:pointer;transition:all 0.2s;'",
      " title='", title_txt, "'",
      " onmouseover=\"this.style.backgroundColor='#f5f5f5'\"",
      " onmouseout=\"this.style.backgroundColor='rgba(255,255,255,0.95)'\"",
      " onclick=\"", onclick_fn, "\">",
      "<svg width='60%' height='60%' viewBox='0 0 24 24' fill='none' stroke='#333'",
      " stroke-width='2.5' stroke-linecap='round' stroke-linejoin='round'>",
      svg_body, "</svg></div>"
    ))
  }

  ui <- htmltools::tagList(
    htmltools::tags$style(htmltools::HTML(css)),
    htmltools::tags$script(htmltools::HTML(js)),
    
    htmltools::tags$div(class = "wsim-container",
      
      # LEFT SIDEBAR
      htmltools::tags$div(id = "wsim_sidebar"),
      
      # RIGHT MAIN AREA
      htmltools::tags$div(class = "wsim-main",
        
        # Tabs Header
        htmltools::tags$div(class = "wsim-tab-header",
          htmltools::tags$button(
            class = "wsim-tab-btn wsim-active",
            onclick = "openWsimTab(event, 'wsim_tab_network')",
            t$network_view
          ),
          htmltools::tags$button(
            class = "wsim-tab-btn",
            onclick = "openWsimTab(event, 'wsim_tab_pcsd')",
            t$pcsd_chart
          )
        ),
        
        # Tab 1 - Network
        htmltools::tags$div(id = "wsim_tab_network", class = "wsim-tab-content wsim-active-tab",
          style = "display:flex; position:relative;",
          network_plot
        ),
        
        # Tab 2 - PCSD (Plotly)
        htmltools::tags$div(id = "wsim_tab_pcsd", class = "wsim-tab-content",
          style = "position:relative;",
          htmltools::tags$div(id = "wsim_pcsd_plot", style = "flex:1; width:100%; height:100%; position:relative;", plotly_plot),
          htmltools::HTML("<div style='position:absolute;bottom:15px;right:15px;z-index:1000;display:flex;flex-direction:column;gap:8px;align-items:center;'>"),
          .btn(0, t$download_img, "downloadSimPcsd()", "<path d='M21 15v4a2 2 0 0 1-2 2H5a2 2 0 0 1-2-2v-4'></path><polyline points='7 10 12 15 17 10'></polyline><line x1='12' y1='15' x2='12' y2='3'></line>"),
          .btn(0, t$info, "var m=document.getElementById('wsim_info_pcsd'); m.style.display=(m.style.display==='block'?'none':'block');", "<circle cx='12' cy='12' r='10'></circle><line x1='12' y1='16' x2='12' y2='12'></line><line x1='12' y1='8' x2='12.01' y2='8'></line>"),
          .btn(0, t$fullscreen, fs_onclick, svg_fs),
          htmltools::HTML("</div>"),
          
          # PCSD Info Modal
          htmltools::HTML(paste0("
          <div id='wsim_info_pcsd' style='position: absolute; top: 50%; left: 50%; transform: translate(-50%, -50%); width: 90%; max-width: 450px; max-height: 80vh; overflow-y: auto; box-sizing: border-box; background-color: #fff; z-index: 2000; padding: 20px; border-radius: 8px; box-shadow: 0 4px 20px rgba(0,0,0,0.2); border: 1px solid #eaeaea; display: none; font-family: Inter, Roboto, sans-serif;'>
            <div style='display:flex; justify-content:space-between; align-items:center; border-bottom:1px solid #eaeaea; padding-bottom:10px; margin-bottom:15px;'>
              <h3 style='margin:0; color:#444; font-size:16px;'>", t$pcsd_chart, "</h3>
              <span style='cursor:pointer; font-size:20px; font-weight:bold; color:#888; line-height:1;' onclick=\"document.getElementById('wsim_info_pcsd').style.display='none';\">&times;</span>
            </div>
            <div style='margin:0; color:#666; font-size:13px; line-height:1.6;'>", t$info_text_sim_pcsd, "</div>
          </div>
          "))
        )
      )
    )
  )

  return(htmltools::browsable(ui))
}


# widget_implications -----------------------------------------------------------

#' Implications Widget for Psychlab
#'
#' @description Creates a three-tab HTML widget: Tab 1 shows the Ideal Digraph,
#'   Tab 2 shows the Impact and Feedback Barplot, and Tab 3 shows the Hypothetical
#'   Scenarios Plot.
#'
#' @param x A \code{wimp} object.
#' @param lang Language for the UI. \code{"en"} (default) or \code{"es"}.
#' @param ... Additional arguments passed to the plotting functions.
#'
#' @return A \code{browsable} HTML object with a three-tab interface.
#' @export
#'
#' @importFrom htmltools tagList tags browsable HTML
#'
#' @examples
#' \dontrun{
#' widget_implications(example_wimp, lang = "es")
#' }
widget_implications <- function(x, lang = "en", ...) {
  if (!inherits(x, "wimp")) stop("Input must be a 'wimp' object.")
  if (!lang %in% c("en", "es")) lang <- "en"

  # Translations for the tabs (fallback inline)
  tab_ideal <- ifelse(lang == "es", "Grafo del Ideal", "Ideal Digraph")
  tab_if <- ifelse(lang == "es", "Impacto y Feedback", "Impact & Feedback")
  tab_hypo <- ifelse(lang == "es", "Escenarios Hip.", "Hypo. Scenarios")
  fullscreen_title <- ifelse(lang == "es", "Pantalla Completa", "Fullscreen")
  
  # 1. Generate individual plots
  # Ideal Digraph
  plot_ideal <- idealdigraph(x, lang = lang, height = "100%", width = "100%", ...)
  plot_ideal$width <- "100%"
  plot_ideal$height <- "100%"
  plot_ideal$sizingPolicy$defaultWidth <- "100%"
  plot_ideal$sizingPolicy$defaultHeight <- "100%"
  
  # IF Barchart
  plot_if <- if_barchart(x, lang = lang, ...)
  plot_if$width <- "100%"
  plot_if$height <- "100%"
  plot_if$sizingPolicy$defaultWidth <- "100%"
  plot_if$sizingPolicy$defaultHeight <- "100%"
  
  # Hypo Plot
  plot_hypo <- hypo_plot(x, lang = lang, ...)
  plot_hypo$width <- "100%"
  plot_hypo$height <- "100%"
  plot_hypo$sizingPolicy$defaultWidth <- "100%"
  plot_hypo$sizingPolicy$defaultHeight <- "100%"

  # 2. Custom CSS and JS for 3-tab layout
  css <- "
    body, html { margin: 0; padding: 0; width: 100%; height: 100%; overflow: hidden; }
    * { box-sizing: border-box; }
    .wt-tab-container {
      width: 100%; height: 100%; display: flex; flex-direction: column;
      font-family: 'Inter', 'Roboto', 'Segoe UI', sans-serif;
    }
    .wt-tab-header {
      overflow: hidden; border: 1px solid #eaeaea; background-color: #ffffff;
      display: flex; border-radius: 4px 4px 0 0;
    }
    .wt-tab-header button {
      background-color: inherit; border: none; outline: none; cursor: pointer;
      padding: 12px 16px; transition: 0.3s; font-size: 14px; font-weight: 600;
      color: #666666; flex-grow: 1; text-align: center;
    }
    .wt-tab-header button:hover { background-color: #f9f9f9; }
    .wt-tab-header button.active {
      background-color: #ffffff; color: #8cc63f; border-bottom: 3px solid #8cc63f;
    }
    .wt-tab-content {
      display: none; padding: 0; border: 1px solid #eaeaea; border-top: none;
      flex-grow: 1; height: 100%; width: 100%; position: relative; overflow: hidden;
      border-radius: 0 0 4px 4px;
    }
    .wt-tab-content.active { display: flex; flex-direction: column; }
  "
  css <- paste0(css, .wt_responsive_css())

  js <- "
    function openImplicationsTab(evt, tabName) {
      var i, tabcontent, tablinks;
      
      tabcontent = document.getElementsByClassName('wt-tab-content');
      for (i = 0; i < tabcontent.length; i++) {
        tabcontent[i].style.display = 'none';
        tabcontent[i].classList.remove('active');
      }
      
      tablinks = document.getElementsByClassName('wt-tab-btn');
      for (i = 0; i < tablinks.length; i++) {
        tablinks[i].className = tablinks[i].className.replace(' active', '');
      }
      
      var selectedTab = document.getElementById(tabName);
      selectedTab.style.display = 'flex';
      selectedTab.classList.add('active');
      evt.currentTarget.className += ' active';
      
      // Trigger resize for htmlwidgets
      window.dispatchEvent(new Event('resize'));
    }

    function downloadImplPlot(containerId, filename) {
      var plotEl = document.querySelector('#' + containerId + ' .js-plotly-plot');
      if (plotEl && window.Plotly) {
        Plotly.downloadImage(plotEl, {format: 'png', width: plotEl.clientWidth, height: plotEl.clientHeight, filename: filename});
      }
    }
  "

  # Button Helper
  .btn <- function(title_txt, onclick_fn, svg_body) {
    htmltools::HTML(paste0(
      "<div style='background-color:rgba(255,255,255,0.95);width:clamp(26px, 4vmin, 34px);height:clamp(26px, 4vmin, 34px);border-radius:6px;",
      "box-shadow:0 2px 10px rgba(0,0,0,0.1);border:1px solid #ddd;display:flex;",
      "align-items:center;justify-content:center;cursor:pointer;transition:all 0.2s;'",
      " title='", title_txt, "'",
      " onmouseover=\"this.style.backgroundColor='#f5f5f5'\"",
      " onmouseout=\"this.style.backgroundColor='rgba(255,255,255,0.95)'\"",
      " onclick=\"", onclick_fn, "\">",
      "<svg width='60%' height='60%' viewBox='0 0 24 24' fill='none' stroke='#333'",
      " stroke-width='2.5' stroke-linecap='round' stroke-linejoin='round'>",
      svg_body, "</svg></div>"
    ))
  }

  svg_fs <- "<path d='M8 3H5a2 2 0 0 0-2 2v3m18 0V5a2 2 0 0 0-2-2h-3m0 18h3a2 2 0 0 0 2-2v-3M3 16v3a2 2 0 0 0 2 2h3'></path>"
  svg_info <- "<circle cx='12' cy='12' r='10'></circle><line x1='12' y1='16' x2='12' y2='12'></line><line x1='12' y1='8' x2='12.01' y2='8'></line>"
  svg_down <- "<path d='M21 15v4a2 2 0 0 1-2 2H5a2 2 0 0 1-2-2v-4'></path><polyline points='7 10 12 15 17 10'></polyline><line x1='12' y1='15' x2='12' y2='3'></line>"
  svg_set  <- "<circle cx='12' cy='12' r='3'></circle><path d='M19.4 15a1.65 1.65 0 0 0 .33 1.82l.06.06a2 2 0 0 1 0 2.83 2 2 0 0 1-2.83 0l-.06-.06a1.65 1.65 0 0 0-1.82-.33 1.65 1.65 0 0 0-1 1.51V21a2 2 0 0 1-2 2 2 2 0 0 1-2-2v-.09A1.65 1.65 0 0 0 9 19.4a1.65 1.65 0 0 0-1.82.33l-.06.06a2 2 0 0 1-2.83 0 2 2 0 0 1 0-2.83l.06-.06a1.65 1.65 0 0 0 .33-1.82 1.65 1.65 0 0 0-1.51-1H3a2 2 0 0 1-2-2 2 2 0 0 1 2-2h.09A1.65 1.65 0 0 0 4.6 9a1.65 1.65 0 0 0-.33-1.82l-.06-.06a2 2 0 0 1 0-2.83 2 2 0 0 1 2.83 0l.06.06a1.65 1.65 0 0 0 1.82.33H9a1.65 1.65 0 0 0 1-1.51V3a2 2 0 0 1 2-2 2 2 0 0 1 2 2v.09a1.65 1.65 0 0 0 1 1.51 1.65 1.65 0 0 0 1.82-.33l.06-.06a2 2 0 0 1 2.83 0 2 2 0 0 1 0 2.83l-.06.06a1.65 1.65 0 0 0-.33 1.82V9a1.65 1.65 0 0 0 1.51 1H21a2 2 0 0 1 2 2 2 2 0 0 1-2 2h-.09a1.65 1.65 0 0 0-1.51 1z'></path>"

  t <- wt_i18n(lang)

  fs_onclick <- "var el=this.closest('.wt-tab-container')||this.closest('.wt-tab-content'); if(!document.fullscreenElement){el.requestFullscreen().catch(e=>console.log(e))}else{document.exitFullscreen()}"
  
  gen_buttons <- function(tab_id) {
    htmltools::HTML(paste0(
      "<div style='position:absolute;bottom:15px;right:15px;z-index:1000;display:flex;flex-direction:column;gap:8px;align-items:center;'>",
      .btn(t$info, paste0("var m=document.getElementById('", tab_id, "_info_modal'); m.style.display=(m.style.display==='block'?'none':'block');"), svg_info),
      .btn(t$hm_settings, paste0("var m=document.getElementById('", tab_id, "_settings_modal'); m.style.display=(m.style.display==='block'?'none':'block');"), svg_set),
      .btn(t$export_png, paste0("downloadImplPlot('", tab_id, "', '", tab_id, "_Export');"), svg_down),
      .btn(fullscreen_title, fs_onclick, svg_fs),
      "</div>"
    ))
  }
  
  gen_modals <- function(tab_id, title) {
    htmltools::HTML(paste0("
      <div id='", tab_id, "_info_modal' style='position: absolute; top: 50%; left: 50%; transform: translate(-50%, -50%); width: 90%; max-width: 450px; max-height: 80vh; overflow-y: auto; box-sizing: border-box; background-color: #fff; z-index: 2000; padding: 20px; border-radius: 8px; box-shadow: 0 4px 20px rgba(0,0,0,0.2); border: 1px solid #eaeaea; display: none; font-family: Inter, Roboto, sans-serif;'>
        <div style='display:flex; justify-content:space-between; align-items:center; border-bottom:1px solid #eaeaea; padding-bottom:10px; margin-bottom:15px;'>
          <h3 style='margin:0; color:#444; font-size:16px;'>", title, " Info</h3>
          <span style='cursor:pointer; font-size:20px; font-weight:bold; color:#888; line-height:1;' onclick=\"document.getElementById('", tab_id, "_info_modal').style.display='none';\">&times;</span>
        </div>
        <p style='margin:0; color:#666; font-size:13px; line-height:1.6;'>Información sobre ", title, "</p>
      </div>
      <div id='", tab_id, "_settings_modal' style='position: absolute; top: 50%; left: 50%; transform: translate(-50%, -50%); width: 90%; max-width: 450px; max-height: 80vh; overflow-y: auto; box-sizing: border-box; background-color: #fff; z-index: 2000; padding: 20px; border-radius: 8px; box-shadow: 0 4px 20px rgba(0,0,0,0.2); border: 1px solid #eaeaea; display: none; font-family: Inter, Roboto, sans-serif;'>
        <div style='display:flex; justify-content:space-between; align-items:center; border-bottom:1px solid #eaeaea; padding-bottom:10px; margin-bottom:15px;'>
          <h3 style='margin:0; color:#444; font-size:16px;'>Ajustes de ", title, "</h3>
          <span style='cursor:pointer; font-size:20px; font-weight:bold; color:#888; line-height:1;' onclick=\"document.getElementById('", tab_id, "_settings_modal').style.display='none';\">&times;</span>
        </div>
        <p style='margin:0; color:#666; font-size:13px; line-height:1.6;'>No hay ajustes configurables en esta versión.</p>
      </div>
    "))
  }

  # 4. Construct UI
  ui <- htmltools::tagList(
    htmltools::tags$style(htmltools::HTML(css)),
    htmltools::tags$script(htmltools::HTML(js)),
    
    htmltools::tags$div(class = "wt-tab-container",
      htmltools::tags$div(class = "wt-tab-header",
        htmltools::tags$button(
          class = "wt-tab-btn active", 
          onclick = "openImplicationsTab(event, 'impl_tab_ideal')", 
          tab_ideal
        ),
        htmltools::tags$button(
          class = "wt-tab-btn", 
          onclick = "openImplicationsTab(event, 'impl_tab_if')", 
          tab_if
        ),
        htmltools::tags$button(
          class = "wt-tab-btn", 
          onclick = "openImplicationsTab(event, 'impl_tab_hypo')", 
          tab_hypo
        )
      ),
      
      # Tab 1: Ideal Digraph
      htmltools::tags$div(id = "impl_tab_ideal", class = "wt-tab-content active", 
        htmltools::tags$div(style = "flex: 1; width: 100%; height: 100%; min-height: 0; position: relative;", plot_ideal)
        # idealdigraph already injects its own floating buttons!
      ),
      
      # Tab 2: IF Barchart
      htmltools::tags$div(id = "impl_tab_if", class = "wt-tab-content", 
        htmltools::tags$div(style = "flex: 1; width: 100%; height: 100%; min-height: 0; position: relative;", plot_if),
        gen_buttons("impl_tab_if"),
        gen_modals("impl_tab_if", tab_if)
      ),
      
      # Tab 3: Hypo Plot
      htmltools::tags$div(id = "impl_tab_hypo", class = "wt-tab-content", 
        htmltools::tags$div(style = "flex: 1; width: 100%; height: 100%; min-height: 0; position: relative;", plot_hypo),
        gen_buttons("impl_tab_hypo"),
        gen_modals("impl_tab_hypo", tab_hypo)
      )
    )
  )
  
  return(htmltools::browsable(ui))
}

#' Generate Tabbed HTML Widget for Wellness Analysis (Análisis de Bienestar) / Adjustment
#'
#' @description
#' Creates a standalone HTML widget containing a two-tab interface for analyzing wellness.
#' If one WimpGrid is provided, it shows the Self analysis and SSI structure.
#' If two WimpGrids are provided, it shows the Self monitoring and SSI monitoring.
#'
#' @param x A `wimp` object representing the baseline evaluation.
#' @param y An optional `wimp` object representing the post evaluation.
#' @param lang Language parameter ("en" or "es").
#' @param ... Additional arguments passed to the underlying plot functions.
#'
#' @return An object of class `htmltools::browsable`.
#' @export
#'
#' @examples
#' \dontrun{
#' widget_adjustment(example_wimp)
#' widget_adjustment(example_wimp, example_wimp_post)
#' }
widget_adjustment <- function(x, y = NULL, lang = "en", ...) {
  if (!inherits(x, "wimp")) stop("Input x must be a 'wimp' object.")
  if (!is.null(y) && !inherits(y, "wimp")) stop("Input y must be a 'wimp' object.")
  if (!lang %in% c("en", "es")) lang <- "en"
  
  t <- wt_i18n(lang)

  tab1_id <- "bienestar_tab1"
  tab2_id <- "bienestar_tab2"

  if (is.null(y)) {
    plot1 <- self_plot(x, ...)
    plot2 <- ssi_heatmap(x, ...)
    tab1_title <- if(lang=="es") "Análisis del Self" else "Self Analysis"
    tab2_title <- if(lang=="es") "Estructura SSI" else "SSI Structure"
  } else {
    plot1 <- monitoring_self(x, y, ...)
    plot2 <- monitoring_ssi(x, y, ...)
    tab1_title <- if(lang=="es") "Monitorización del Self" else "Self Monitoring"
    tab2_title <- if(lang=="es") "Monitorización SSI" else "SSI Monitoring"
  }
  
  # Ensure plots occupy full container (plotly standard)
  if (inherits(plot1, "plotly")) {
    plot1 <- plotly::layout(plot1, autosize = TRUE)
    plot1$sizingPolicy$defaultHeight <- "100%"
    plot1$sizingPolicy$defaultWidth <- "100%"
    plot1$height <- "100%"
    plot1$width <- "100%"
  }
  if (inherits(plot2, "plotly")) {
    plot2 <- plotly::layout(plot2, autosize = TRUE)
    plot2$sizingPolicy$defaultHeight <- "100%"
    plot2$sizingPolicy$defaultWidth <- "100%"
    plot2$height <- "100%"
    plot2$width <- "100%"
  }
  
  css <- "
    .wt-tab-container { width: 100%; height: 100vh; display: flex; flex-direction: column; font-family: 'Inter', Roboto, sans-serif; background: #fafafa; }
    .wt-tab-header { display: flex; background: #fff; border-bottom: 2px solid #eaeaea; padding: 0 10px; flex-shrink: 0; }
    .wt-tab-btn { background: none; border: none; padding: 14px 20px; cursor: pointer; font-size: 14px; font-weight: 600; color: #888; border-bottom: 3px solid transparent; transition: all 0.2s; }
    .wt-tab-btn:hover { color: #333; }
    .wt-tab-btn.active { color: #8cc63f; border-bottom-color: #8cc63f; }
    .wt-tab-content { display: none; flex: 1; min-height: 0; position: relative; }
    .wt-tab-content.active { display: flex; }
  "
  css <- paste0(css, .wt_responsive_css())

  js <- "
    function openBienestarTab(evt, tabId) {
      var i, tabcontent, tablinks;
      tabcontent = document.getElementsByClassName('wt-tab-content');
      for (i = 0; i < tabcontent.length; i++) {
        tabcontent[i].style.display = 'none';
        tabcontent[i].classList.remove('active');
      }
      
      tablinks = document.getElementsByClassName('wt-tab-btn');
      for (i = 0; i < tablinks.length; i++) {
        tablinks[i].className = tablinks[i].className.replace(' active', '');
      }
      
      var selectedTab = document.getElementById(tabId);
      selectedTab.style.display = 'flex';
      selectedTab.classList.add('active');
      evt.currentTarget.className += ' active';
      
      window.dispatchEvent(new Event('resize'));
    }

    function downloadBienestarPlot(containerId, filename) {
      var plotEl = document.querySelector('#' + containerId + ' .js-plotly-plot');
      if (plotEl && window.Plotly) {
        Plotly.downloadImage(plotEl, {format: 'png', width: plotEl.clientWidth, height: plotEl.clientHeight, filename: filename});
      }
    }
  "

  # Button Helper
  .btn <- function(title_txt, onclick_fn, svg_body) {
    htmltools::HTML(paste0(
      "<div style='background-color:rgba(255,255,255,0.95);width:clamp(26px, 4vmin, 34px);height:clamp(26px, 4vmin, 34px);border-radius:6px;",
      "box-shadow:0 2px 10px rgba(0,0,0,0.1);border:1px solid #ddd;display:flex;",
      "align-items:center;justify-content:center;cursor:pointer;transition:all 0.2s;'",
      " title='", title_txt, "'",
      " onmouseover=\"this.style.backgroundColor='#f5f5f5'\"",
      " onmouseout=\"this.style.backgroundColor='rgba(255,255,255,0.95)'\"",
      " onclick=\"", onclick_fn, "\">",
      "<svg width='60%' height='60%' viewBox='0 0 24 24' fill='none' stroke='#333'",
      " stroke-width='2.5' stroke-linecap='round' stroke-linejoin='round'>",
      svg_body, "</svg></div>"
    ))
  }

  svg_fs <- "<path d='M8 3H5a2 2 0 0 0-2 2v3m18 0V5a2 2 0 0 0-2-2h-3m0 18h3a2 2 0 0 0 2-2v-3M3 16v3a2 2 0 0 0 2 2h3'></path>"
  svg_info <- "<circle cx='12' cy='12' r='10'></circle><line x1='12' y1='16' x2='12' y2='12'></line><line x1='12' y1='8' x2='12.01' y2='8'></line>"
  svg_down <- "<path d='M21 15v4a2 2 0 0 1-2 2H5a2 2 0 0 1-2-2v-4'></path><polyline points='7 10 12 15 17 10'></polyline><line x1='12' y1='15' x2='12' y2='3'></line>"
  svg_set  <- "<circle cx='12' cy='12' r='3'></circle><path d='M19.4 15a1.65 1.65 0 0 0 .33 1.82l.06.06a2 2 0 0 1 0 2.83 2 2 0 0 1-2.83 0l-.06-.06a1.65 1.65 0 0 0-1.82-.33 1.65 1.65 0 0 0-1 1.51V21a2 2 0 0 1-2 2 2 2 0 0 1-2-2v-.09A1.65 1.65 0 0 0 9 19.4a1.65 1.65 0 0 0-1.82.33l-.06.06a2 2 0 0 1-2.83 0 2 2 0 0 1 0-2.83l.06.06a1.65 1.65 0 0 0 .33-1.82 1.65 1.65 0 0 0-1.51-1H3a2 2 0 0 1-2-2 2 2 0 0 1 2-2h.09A1.65 1.65 0 0 0 4.6 9a1.65 1.65 0 0 0-.33-1.82l-.06-.06a2 2 0 0 1 0-2.83 2 2 0 0 1 2.83 0l.06.06a1.65 1.65 0 0 0 1.82.33H9a1.65 1.65 0 0 0 1-1.51V3a2 2 0 0 1 2-2 2 2 0 0 1 2 2v.09a1.65 1.65 0 0 0 1 1.51 1.65 1.65 0 0 0 1.82-.33l.06-.06a2 2 0 0 1 2.83 0 2 2 0 0 1 0 2.83l-.06.06a1.65 1.65 0 0 0-.33 1.82V9a1.65 1.65 0 0 0 1.51 1H21a2 2 0 0 1 2 2 2 2 0 0 1-2 2h-.09a1.65 1.65 0 0 0-1.51 1z'></path>"

  fullscreen_title <- if (lang == "es") "Pantalla Completa" else "Fullscreen"
  fs_onclick <- "var el=this.closest('.wt-tab-container')||this.closest('.wt-tab-content'); if(!document.fullscreenElement){el.requestFullscreen().catch(e=>console.log(e))}else{document.exitFullscreen()}"
  
  gen_buttons <- function(tab_id) {
    htmltools::HTML(paste0(
      "<div style='position:absolute;bottom:15px;right:15px;z-index:1000;display:flex;flex-direction:column;gap:8px;align-items:center;'>",
      .btn(t$info, paste0("var m=document.getElementById('", tab_id, "_info_modal'); m.style.display=(m.style.display==='block'?'none':'block');"), svg_info),
      .btn(t$hm_settings, paste0("var m=document.getElementById('", tab_id, "_settings_modal'); m.style.display=(m.style.display==='block'?'none':'block');"), svg_set),
      .btn(t$export_png, paste0("downloadBienestarPlot('", tab_id, "', '", tab_id, "_Export');"), svg_down),
      .btn(fullscreen_title, fs_onclick, svg_fs),
      "</div>"
    ))
  }
  
  gen_modals <- function(tab_id, title) {
    htmltools::HTML(paste0("
      <div id='", tab_id, "_info_modal' style='position: absolute; top: 50%; left: 50%; transform: translate(-50%, -50%); width: 90%; max-width: 450px; max-height: 80vh; overflow-y: auto; box-sizing: border-box; background-color: #fff; z-index: 2000; padding: 20px; border-radius: 8px; box-shadow: 0 4px 20px rgba(0,0,0,0.2); border: 1px solid #eaeaea; display: none; font-family: Inter, Roboto, sans-serif;'>
        <div style='display:flex; justify-content:space-between; align-items:center; border-bottom:1px solid #eaeaea; padding-bottom:10px; margin-bottom:15px;'>
          <h3 style='margin:0; color:#444; font-size:16px;'>", title, " Info</h3>
          <span style='cursor:pointer; font-size:20px; font-weight:bold; color:#888; line-height:1;' onclick=\"document.getElementById('", tab_id, "_info_modal').style.display='none';\">&times;</span>
        </div>
        <p style='margin:0; color:#666; font-size:13px; line-height:1.6;'>Información sobre ", title, "</p>
      </div>
      <div id='", tab_id, "_settings_modal' style='position: absolute; top: 50%; left: 50%; transform: translate(-50%, -50%); width: 90%; max-width: 450px; max-height: 80vh; overflow-y: auto; box-sizing: border-box; background-color: #fff; z-index: 2000; padding: 20px; border-radius: 8px; box-shadow: 0 4px 20px rgba(0,0,0,0.2); border: 1px solid #eaeaea; display: none; font-family: Inter, Roboto, sans-serif;'>
        <div style='display:flex; justify-content:space-between; align-items:center; border-bottom:1px solid #eaeaea; padding-bottom:10px; margin-bottom:15px;'>
          <h3 style='margin:0; color:#444; font-size:16px;'>Ajustes de ", title, "</h3>
          <span style='cursor:pointer; font-size:20px; font-weight:bold; color:#888; line-height:1;' onclick=\"document.getElementById('", tab_id, "_settings_modal').style.display='none';\">&times;</span>
        </div>
        <p style='margin:0; color:#666; font-size:13px; line-height:1.6;'>No hay ajustes configurables en esta versión.</p>
      </div>
    "))
  }

  ui <- htmltools::tagList(
    htmltools::tags$style(htmltools::HTML(css)),
    htmltools::tags$script(htmltools::HTML(js)),
    
    htmltools::tags$div(class = "wt-tab-container",
      htmltools::tags$div(class = "wt-tab-header",
        htmltools::tags$button(
          class = "wt-tab-btn active", 
          onclick = "openBienestarTab(event, 'bienestar_tab1')", 
          tab1_title
        ),
        htmltools::tags$button(
          class = "wt-tab-btn", 
          onclick = "openBienestarTab(event, 'bienestar_tab2')", 
          tab2_title
        )
      ),
      
      # Tab 1
      htmltools::tags$div(id = "bienestar_tab1", class = "wt-tab-content active", 
        htmltools::tags$div(style = "flex: 1; width: 100%; height: 100%; min-height: 0; position: relative;", plot1),
        gen_buttons("bienestar_tab1"),
        gen_modals("bienestar_tab1", tab1_title)
      ),
      
      # Tab 2
      htmltools::tags$div(id = "bienestar_tab2", class = "wt-tab-content", 
        htmltools::tags$div(style = "flex: 1; width: 100%; height: 100%; min-height: 0; position: relative;", plot2),
        gen_buttons("bienestar_tab2"),
        gen_modals("bienestar_tab2", tab2_title)
      )
    )
  )
  
  return(htmltools::browsable(ui))
}

# widget_repgrid_biplot -------------------------------------------------------

#' RepGrid Biplot Widget for Psychlab
#'
#' @description Creates a two-tab HTML widget: Tab 1 shows the 2D biplot of a
#'   repertory grid, Tab 2 shows the 3D biplot.
#'
#' @param x An \code{OpenRepGrid} \code{repgrid} object or a numeric ratings
#'   matrix (see \code{\link{repgrid_biplot}}).
#' @param lang Language for the UI. \code{"en"} (default) or \code{"es"}.
#' @param ... Additional arguments passed to \code{repgrid_biplot}.
#'
#' @return A \code{browsable} HTML object with a two-tab interface.
#' @export
#'
#' @importFrom htmltools tagList tags browsable HTML
widget_repgrid_biplot <- function(x, lang = "en", ...) {
  if (!lang %in% c("en", "es")) lang <- "en"
  t <- wt_i18n(lang)
  plots <- list(
    repgrid_biplot(x, dim = 2, ...) %>% plotly::config(displayModeBar = FALSE),
    repgrid_biplot(x, dim = 3, ...) %>% plotly::config(displayModeBar = FALSE)
  )
  .rg_tabbed_widget(plots, c(t$biplot_2d_tab, t$biplot_3d_tab),
                    c("RepGrid_Biplot_2D", "RepGrid_Biplot_3D"),
                    t$info_text_biplot, t)
}

# Internal: tabbed plotly widget shared by the RepGrid widgets ----------------
.rg_tabbed_widget <- function(plots, labels, fnames, info, t, settings = NULL,
                              export_js = NULL, extra_js = NULL) {
  css <- "
    body, html { margin: 0; padding: 0; width: 100%; height: 100%; overflow: hidden; }
    * { box-sizing: border-box; }
    .wrg-container { position: relative; width: 100%; height: 100%; display: flex; flex-direction: column;
      font-family: 'Inter', 'Roboto', 'Segoe UI', sans-serif; }
    .wrg-tab-header { overflow: hidden; border: 1px solid #eaeaea; background-color: #ffffff;
      display: flex; border-radius: 4px 4px 0 0; }
    .wrg-tab-header button { background-color: inherit; border: none; outline: none; cursor: pointer;
      padding: 12px 24px; transition: 0.3s; font-size: 14px; font-weight: 600; color: #666666; flex-grow: 1; }
    .wrg-tab-header button:hover { background-color: #f9f9f9; }
    .wrg-tab-header button.wrg-active { background-color: #ffffff; color: #8cc63f; border-bottom: 2px solid #8cc63f; }
    .wrg-tab-content { position: relative; display: none; padding: 0; border: 1px solid #ccc;
      border-top: none; flex-grow: 1; height: 0; background-color: #ffffff; border-radius: 0 0 4px 4px; }
    .wrg-tab-content > .html-widget { width: 100% !important; height: 100% !important; flex-grow: 1; }
  "
  css <- paste0(css, .wt_responsive_css())

  js <- "
    function openWrgTab(evt, tabName) {
      var i, c = document.getElementsByClassName('wrg-tab-content'), b = document.getElementsByClassName('wrg-tab-btn');
      for (i = 0; i < c.length; i++) { c[i].style.display = 'none'; }
      for (i = 0; i < b.length; i++) { b[i].className = b[i].className.replace(' wrg-active', ''); }
      document.getElementById(tabName).style.display = 'flex';
      document.getElementById(tabName).style.flexDirection = 'column';
      evt.currentTarget.className += ' wrg-active';
      window.dispatchEvent(new Event('resize'));
    }
    function downloadWrgPlot(tabId, name) {
      var el = document.querySelector('#' + tabId + ' .js-plotly-plot');
      if (el && window.Plotly) {
        Plotly.downloadImage(el, {format: 'png', width: el.clientWidth, height: el.clientHeight, filename: name});
      }
    }
  "

  svg_dl <- "<path d='M21 15v4a2 2 0 0 1-2 2H5a2 2 0 0 1-2-2v-4'></path><polyline points='7 10 12 15 17 10'></polyline><line x1='12' y1='15' x2='12' y2='3'></line>"
  svg_i  <- "<circle cx='12' cy='12' r='10'></circle><line x1='12' y1='16' x2='12' y2='12'></line><line x1='12' y1='8' x2='12.01' y2='8'></line>"
  svg_fs <- "<path d='M8 3H5a2 2 0 0 0-2 2v3m18 0V5a2 2 0 0 0-2-2h-3m0 18h3a2 2 0 0 0 2-2v-3M3 16v3a2 2 0 0 0 2 2h3'></path>"

  svg_set <- "<circle cx='12' cy='12' r='3'></circle><path d='M19.4 15a1.65 1.65 0 0 0 .33 1.82l.06.06a2 2 0 0 1 0 2.83 2 2 0 0 1-2.83 0l-.06-.06a1.65 1.65 0 0 0-1.82-.33 1.65 1.65 0 0 0-1 1.51V21a2 2 0 0 1-2 2 2 2 0 0 1-2-2v-.09A1.65 1.65 0 0 0 9 19.4a1.65 1.65 0 0 0-1.82.33l-.06.06a2 2 0 0 1-2.83 0 2 2 0 0 1 0-2.83l.06-.06a1.65 1.65 0 0 0 .33-1.82 1.65 1.65 0 0 0-1.51-1H3a2 2 0 0 1-2-2 2 2 0 0 1 2-2h.09A1.65 1.65 0 0 0 4.6 9a1.65 1.65 0 0 0-.33-1.82l-.06-.06a2 2 0 0 1 0-2.83 2 2 0 0 1 2.83 0l.06.06a1.65 1.65 0 0 0 1.82.33H9a1.65 1.65 0 0 0 1-1.51V3a2 2 0 0 1 2-2 2 2 0 0 1 2 2v.09a1.65 1.65 0 0 0 1 1.51 1.65 1.65 0 0 0 1.82-.33l.06-.06a2 2 0 0 1 2.83 0 2 2 0 0 1 0 2.83l-.06.06a1.65 1.65 0 0 0-.33 1.82V9a1.65 1.65 0 0 0 1.51 1H21a2 2 0 0 1 2 2 2 2 0 0 1-2 2h-.09a1.65 1.65 0 0 0-1.51 1z'></path>"

  .btn <- function(title_txt, onclick_fn, svg_body) {
    htmltools::HTML(paste0(
      "<div style='background-color:rgba(255,255,255,0.95);width:clamp(26px, 4vmin, 34px);height:clamp(26px, 4vmin, 34px);border-radius:6px;",
      "box-shadow:0 2px 10px rgba(0,0,0,0.1);border:1px solid #ddd;display:flex;",
      "align-items:center;justify-content:center;cursor:pointer;transition:all 0.2s;'",
      " title='", title_txt, "'",
      " onmouseover=\"this.style.backgroundColor='#f5f5f5'\"",
      " onmouseout=\"this.style.backgroundColor='rgba(255,255,255,0.95)'\"",
      " onclick=\"", onclick_fn, "\">",
      "<svg width='60%' height='60%' viewBox='0 0 24 24' fill='none' stroke='#333'",
      " stroke-width='2.5' stroke-linecap='round' stroke-linejoin='round'>",
      svg_body, "</svg></div>"))
  }

  .modal <- function(id, title_txt, body_txt) {
    htmltools::HTML(paste0(
      "<div id='", id, "' style='position:absolute;top:50%;left:50%;",
      "transform:translate(-50%,-50%);width:90%;max-width:450px;max-height:80vh;overflow-y:auto;box-sizing:border-box;background:#fff;",
      "z-index:2000;padding:20px;border-radius:8px;",
      "box-shadow:0 4px 20px rgba(0,0,0,0.2);border:1px solid #eaeaea;display:none;",
      "font-family:Inter,Roboto,sans-serif;'>",
      "<div style='display:flex;justify-content:space-between;align-items:center;",
      "border-bottom:1px solid #eaeaea;padding-bottom:10px;margin-bottom:15px;'>",
      "<h3 style='margin:0;color:#444;font-size:16px;'>", title_txt, "</h3>",
      "<span style='cursor:pointer;font-size:20px;font-weight:bold;color:#888;line-height:1;'",
      " onclick=\"document.getElementById('", id, "').style.display='none';\">&times;</span>",
      "</div>",
      "<p style='margin:0;color:#666;font-size:13px;line-height:1.6;'>", body_txt, "</p>",
      "</div>"))
  }

  fs_onclick <- "var el=this.closest('.wrg-container'); if(!document.fullscreenElement){el.requestFullscreen().catch(e=>console.log(e))}else{document.exitFullscreen()}"

  .tab <- function(id, plot, modal_id, tab_title, fname, display, i) {
    htmltools::tags$div(id = id, class = "wrg-tab-content",
      style = if (display) "display:flex; flex-direction:column;" else NULL,
      plot,
      .modal(modal_id, tab_title, info),
      htmltools::HTML("<div style='position:absolute;bottom:15px;right:15px;z-index:1000;display:flex;flex-direction:column;gap:8px;align-items:center;'>"),
      .btn(t$info, sprintf("var m=document.getElementById('%s'); m.style.display=(m.style.display==='block'?'none':'block');", modal_id), svg_i),
      if (!is.null(settings)) .btn(t$hm_settings, "var m=document.getElementById('wrg_settings_modal'); m.style.display=(m.style.display==='block'?'none':'block');", svg_set),
      .btn(if (is.null(export_js) || is.na(export_js[i])) t$export_png else t$export_csv,
           if (is.null(export_js) || is.na(export_js[i])) sprintf("downloadWrgPlot('%s','%s')", id, fname) else export_js[i],
           svg_dl),
      .btn(t$fullscreen, fs_onclick, svg_fs),
      htmltools::HTML("</div>"))
  }

  n <- length(plots)
  ids <- paste0("wrg_tab_", seq_len(n))
  ui <- htmltools::tagList(
    htmltools::tags$style(htmltools::HTML(css)),
    htmltools::tags$script(htmltools::HTML(js)),
    if (!is.null(settings)) htmltools::tags$script(htmltools::HTML(settings$js)),
    if (!is.null(extra_js)) htmltools::tags$script(htmltools::HTML(extra_js)),
    htmltools::tags$div(class = "wrg-container",
      htmltools::tags$div(class = "wrg-tab-header",
        lapply(seq_len(n), function(i) {
          htmltools::tags$button(
            class = if (i == 1) "wrg-tab-btn wrg-active" else "wrg-tab-btn",
            onclick = sprintf("openWrgTab(event, '%s')", ids[i]), labels[i])
        })),
      lapply(seq_len(n), function(i) {
        .tab(ids[i], plots[[i]], paste0("wrg_info_", i), labels[i], fnames[i], i == 1, i)
      }),
      if (!is.null(settings)) htmltools::HTML(paste0(
        "<div id='wrg_settings_modal' style='position:absolute;",
        "top:55px;right:10px;left:auto;transform:none;width:90%;max-width:min(300px, 92vw);max-height:80vh;overflow-y:auto;box-sizing:border-box;background:#fff;",
        "z-index:2000;padding:20px;border-radius:8px;box-shadow:0 4px 20px rgba(0,0,0,0.2);border:1px solid #eaeaea;display:none;",
        "font-family:Inter,Roboto,sans-serif;'>",
        "<div style='display:flex;justify-content:space-between;align-items:center;",
        "border-bottom:1px solid #eaeaea;padding-bottom:10px;margin-bottom:15px;'>",
        "<h3 style='margin:0;color:#444;font-size:16px;'>", t$hm_settings, "</h3>",
        "<span style='cursor:pointer;font-size:20px;font-weight:bold;color:#888;line-height:1;'",
        " onclick=\"document.getElementById('wrg_settings_modal').style.display='none';\">&times;</span>",
        "</div>", settings$html, "</div>"))
    )
  )

  htmltools::browsable(ui)
}

# widget_repgrid_cluster ------------------------------------------------------

#' RepGrid Cluster Widget for Psychlab
#'
#' @description Creates a two-tab HTML widget: Tab 1 shows the dendrogram of
#'   the constructs, Tab 2 the dendrogram of the elements
#'   (see \code{\link{repgrid_cluster}}).
#'
#' @param x An \code{OpenRepGrid} \code{repgrid} object or a numeric ratings
#'   matrix.
#' @param lang Language for the UI. \code{"en"} (default) or \code{"es"}.
#' @param ... Additional arguments passed to \code{repgrid_cluster}.
#'
#' @return A \code{browsable} HTML object with a two-tab interface.
#' @export
widget_repgrid_cluster <- function(x, lang = "en", ...) {
  if (!lang %in% c("en", "es")) lang <- "en"
  t <- wt_i18n(lang)
  args <- list(...)
  dist0   <- if (is.null(args$dist)) "euclidean" else args$dist
  method0 <- if (is.null(args$method)) "ward.D" else args$method
  dists   <- unique(c(.rg_cluster_dists, dist0))
  methods <- unique(c(.rg_cluster_methods, method0))

  plots <- list(
    repgrid_cluster(x, along = "constructs", ...) %>% plotly::config(displayModeBar = FALSE),
    repgrid_cluster(x, along = "elements", ...) %>% plotly::config(displayModeBar = FALSE)
  )

  # Precomputed layouts for every dist x method, swapped in by the menu
  data <- list(
    wrg_tab_1 = .rg_cluster_options(x, "constructs", dists, methods),
    wrg_tab_2 = .rg_cluster_options(x, "elements",   dists, methods)
  )
  data_json <- jsonlite::toJSON(data, auto_unbox = FALSE, na = "null", digits = NA)

  opts <- function(v, sel) paste0(sprintf("<option value='%s'%s>%s</option>", v,
                                  ifelse(v == sel, " selected", ""), v), collapse = "")
  sel_style <- "width:100%;padding:8px 10px;border:1px solid #ddd;border-radius:6px;font-size:13px;font-family:inherit;outline:none;background:#fff;"
  lbl_style <- "display:block;margin:0 0 6px 0;font-size:12px;font-weight:600;color:#555;"
  html <- paste0(
    "<label style='", lbl_style, "'>", t$cluster_dist, "</label>",
    "<select id='wrg_dist' style='", sel_style, "margin-bottom:14px;' onchange='wrgUpdateCluster()'>",
    opts(dists, dist0), "</select>",
    "<label style='", lbl_style, "'>", t$cluster_method, "</label>",
    "<select id='wrg_method' style='", sel_style, "' onchange='wrgUpdateCluster()'>",
    opts(methods, method0), "</select>")

  js <- paste0("
    var WRG_DATA = ", data_json, ";
    function wrgUpdateCluster() {
      var key = document.getElementById('wrg_dist').value + '|' + document.getElementById('wrg_method').value;
      Object.keys(WRG_DATA).forEach(function(tab) {
        var d = WRG_DATA[tab][key], el = document.querySelector('#' + tab + ' .js-plotly-plot');
        if (!d || !el || !window.Plotly) return;
        Plotly.restyle(el, {x: [d.sx], y: [d.sy]}, [0]);
        Plotly.restyle(el, {y: [d.pos]}, [1]);
        Plotly.restyle(el, {x: [d.hx], y: [d.hy], text: [d.hx.map(function(h){ return 'd = ' + h.toFixed(2); })]}, [2]);
        Plotly.relayout(el, {'xaxis.range': [d.hmax * 1.05, -d.hmax * 0.02], 'yaxis.ticktext': d.ticktext});
        // The label set/order just changed (new dist/method): drop the
        // cached full-length labels so wrgFitClusterPlot re-derives them
        // from this new ticktext instead of re-trimming the old one.
        el._wrgOrigTicktext = null;
        wrgFitClusterPlot(el);
      });
    }
  ")

  # Labels here are the y-axis's own ticktext (right side, automargin off),
  # not custom annotations like the dilemmas widget - so instead of
  # repositioning annotation objects, this shrinks the tick font and
  # right margin to the widget's real width, and once that's not enough,
  # truncates the tick strings with an ellipsis (canvas-measured, so it
  # never overflows). The right margin is sized from the CURRENT labels'
  # own longest width rather than a fixed guess, so a tab with short
  # labels (elements) naturally gives the dendrogram far more room than
  # one with long labels (constructs) - the two tabs are no longer forced
  # to share the same oversized margin. The dendrogram's own drawing area
  # (between the two margins) never drops below 200px.
  fit_js <- "
    function wrgClMeasure(text, font) {
      var c = wrgClMeasure._c || (wrgClMeasure._c = document.createElement('canvas'));
      var ctx = c.getContext('2d');
      ctx.font = font;
      return ctx.measureText(text).width;
    }
    function wrgClFitText(text, font, maxWidth) {
      if (maxWidth <= 0) return '';
      if (wrgClMeasure(text, font) <= maxWidth) return text;
      var lo = 0, hi = text.length;
      while (lo < hi) {
        var mid = Math.ceil((lo + hi) / 2);
        var candidate = text.slice(0, mid) + '\\u2026';
        if (wrgClMeasure(candidate, font) <= maxWidth) lo = mid; else hi = mid - 1;
      }
      return lo === 0 ? '' : text.slice(0, lo) + '\\u2026';
    }
    function wrgClPlain(html) { return String(html).replace(/<[^>]*>/g, ''); }

    function wrgFitClusterPlot(el) {
      if (!el || !el.layout || !window.Plotly) return;
      var totalW = document.body.clientWidth; if (!totalW) return;
      if (!el._wrgOrigTicktext) el._wrgOrigTicktext = (el.layout.yaxis.ticktext || []).slice();
      var orig = el._wrgOrigTicktext;

      var iconGutter = totalW <= 340 ? 38 : (totalW <= 520 ? 46 : 0);
      var avail = totalW - iconGutter;
      var marginL = 40;
      var minInner = 200;
      var fontPx = Math.round(Math.max(9, Math.min(12, totalW / 60)));
      var font = 'bold ' + fontPx + 'px \"Open Sans\", verdana, arial, sans-serif';

      var longest = 0;
      orig.forEach(function(txt) {
        var w2 = wrgClMeasure(wrgClPlain(txt), font);
        if (w2 > longest) longest = w2;
      });
      var marginR = Math.min(330, longest + 24);
      if (avail - marginL - marginR < minInner) {
        marginR = Math.max(30, avail - marginL - minInner);
      }
      if (avail - marginL - marginR < minInner) {
        // Margin is already at its floor and it's still not enough: claim
        // back whatever the icon gutter can spare.
        iconGutter = Math.max(0, totalW - marginL - minInner - marginR);
        avail = totalW - iconGutter;
      }
      var elWidth = totalW - iconGutter;
      el.style.setProperty('width', elWidth + 'px', 'important');

      var budget = marginR - 16;
      var newTicktext = orig.map(function(txt) {
        var plain = wrgClPlain(txt);
        if (wrgClMeasure(plain, font) <= budget) return txt;
        // Doesn't fit even at this font size: fall back to the plain,
        // truncated text (losing the bold/colour styling only in that case).
        return wrgClFitText(plain, font, budget);
      });

      Plotly.relayout(el, {
        'yaxis.ticktext': newTicktext,
        'yaxis.tickfont.size': fontPx,
        'margin.l': marginL,
        'margin.r': marginR
      });
    }
    function wrgClusterPlotEls() {
      return [document.querySelector('#wrg_tab_1 .js-plotly-plot'), document.querySelector('#wrg_tab_2 .js-plotly-plot')];
    }
    window.addEventListener('load', function() {
      setTimeout(function() { wrgClusterPlotEls().forEach(wrgFitClusterPlot); }, 50);
      var wrgClLastWidth = document.body.clientWidth;
      setInterval(function() {
        var w = document.body.clientWidth;
        if (w && w !== wrgClLastWidth) {
          wrgClLastWidth = w;
          wrgClusterPlotEls().forEach(wrgFitClusterPlot);
        }
      }, 300);
    });
    window.addEventListener('resize', function() {
      wrgClusterPlotEls().forEach(wrgFitClusterPlot);
    });
  "

  .rg_tabbed_widget(plots, c(t$cluster_constructs_tab, t$cluster_elements_tab),
                    c("RepGrid_Cluster_Constructs", "RepGrid_Cluster_Elements"),
                    t$info_text_cluster, t,
                    settings = list(html = html, js = js),
                    extra_js = fit_js)
}


# widget_repgrid_dilemmas -----------------------------------------------------

#' RepGrid Implicative Dilemmas Widget for Psychlab
#'
#' @description Creates a two-tab HTML widget: Tab 1 shows the diagram of the
#'   implicative dilemmas, Tab 2 a table with each dilemma and the summary
#'   indices (see \code{\link{repgrid_dilemmas}}).
#'
#' @param x An \code{OpenRepGrid} \code{repgrid} object. The self is assumed
#'   in the first column and the ideal in the last one.
#' @param lang Language for the UI. \code{"en"} (default) or \code{"es"}.
#' @param ... Additional arguments passed to \code{OpenRepGrid::indexDilemma}.
#'
#' @return A \code{browsable} HTML object with a two-tab interface.
#' @export
widget_repgrid_dilemmas <- function(x, lang = "en", ...) {
  if (!lang %in% c("en", "es")) lang <- "en"
  t <- wt_i18n(lang)
  args <- list(...)
  mode0 <- if (is.null(args$diff.mode)) 1 else args$diff.mode
  rmin0 <- if (is.null(args$r.min)) 0.35 else args$r.min
  base  <- args[!names(args) %in% c("diff.mode", "r.min")]
  modes <- unique(c(1, 0, mode0))
  rmins <- sort(unique(c(round(seq(0.2, 0.7, by = 0.05), 2), rmin0)))

  esc <- function(v) htmltools::htmlEscape(v)
  rows_html <- function(dl) {
    if (nrow(dl) == 0) {
      paste0("<tr><td colspan='3' style='text-align:center;padding:20px;color:#999;'>", t$no_dilemmas, "</td></tr>")
    } else {
      paste0(sprintf("<tr><td>%s</td><td>%s</td><td style='text-align:right;'>%.3f</td></tr>",
                     esc(dl$congruent), esc(dl$discrepant), dl$r), collapse = "\n")
    }
  }
  view <- function(mode, rmin, inv = FALSE) {
    a <- c(base, list(diff.mode = mode, r.min = rmin))
    dd <- do.call(.rg_dilemma_data, c(list(x), a))
    p  <- do.call(repgrid_dilemmas, c(list(x), a, list(only_involved = inv)))
    b  <- plotly::plotly_build(p)$x
    list(dd = dd, plot = p,
         json = list(data = b$data, layout = b$layout),
         rows = rows_html(dd$dilemmas),
         kpi = c(n = as.character(dd$n_ids), pid = sprintf("%.1f%%", 100 * dd$pid),
                 iid = sprintf("%.1f", dd$iid), picid = sprintf("%.2f", dd$picid)))
  }

  # Precompute every diff.mode x r.min view; the menu swaps them in
  views <- list()
  for (m in modes) for (r in rmins) {
    views[[paste(m, r, "all", sep = "|")]] <- view(m, r, FALSE)
    views[[paste(m, r, "inv", sep = "|")]] <- view(m, r, TRUE)
  }
  cur <- views[[paste(mode0, rmin0, "all", sep = "|")]]

  plot <- cur$plot %>% plotly::config(displayModeBar = FALSE)
  payload <- lapply(views, function(v) list(fig = v$json, rows = v$rows, kpi = as.list(v$kpi)))
  payload_json <- jsonlite::toJSON(payload, auto_unbox = TRUE, null = "null", na = "null", digits = NA)

  kpi <- function(label, id, value) paste0(
    "<div style='flex:1;min-width:110px;border:1px solid #eaeaea;border-radius:8px;padding:10px 14px;background:#fbfdf9;'>",
    "<div style='font-size:11px;color:#888;font-weight:600;'>", label, "</div>",
    "<div id='", id, "' style='font-size:20px;color:#333;font-weight:700;'>", value, "</div></div>")
  table_html <- paste0(
    "<div style='width:100%;height:100%;display:flex;flex-direction:column;font-family:Inter,Roboto,sans-serif;'>",
    "<div style='display:flex;gap:12px;flex-wrap:wrap;padding:16px 20px 8px 20px;flex-shrink:0;'>",
    kpi(t$n_dilemmas, "wrg_kpi_n", cur$kpi[["n"]]), kpi("PID", "wrg_kpi_pid", cur$kpi[["pid"]]),
    kpi("IID", "wrg_kpi_iid", cur$kpi[["iid"]]), kpi("PICID", "wrg_kpi_picid", cur$kpi[["picid"]]),
    "</div>",
    "<div style='flex:1;overflow:auto;padding:8px clamp(16px, 10vw, 70px) 60px 20px;'>",
    "<table id='wrg_dil_table' style='width:100%;border-collapse:collapse;font-size:13px;'>",
    "<thead><tr style='background:#f8f9fa;border-bottom:2px solid #8cc63f;'>",
    "<th>", t$congruent, "</th><th>", t$discrepant, "</th><th style='text-align:right;'>", t$correlation, "</th>",
    "</tr></thead><tbody id='wrg_dil_tbody'>", cur$rows, "</tbody></table></div></div>",
    "<style>#wrg_dil_table thead th{padding:10px 14px;text-align:left;font-weight:600;color:#555;font-size:12px;}",
    "#wrg_dil_table tbody tr{border-bottom:1px solid #f0f0f0;}#wrg_dil_table tbody tr:hover{background:#f4f9ef;}",
    "#wrg_dil_table tbody td{padding:9px 14px;color:#333;}</style>")

  sel_style <- "width:100%;padding:8px 10px;border:1px solid #ddd;border-radius:6px;font-size:13px;font-family:inherit;outline:none;background:#fff;"
  lbl_style <- "display:block;margin:0 0 6px 0;font-size:12px;font-weight:600;color:#555;"
  mode_lab <- c(`1` = t$dil_mode_diff, `0` = t$dil_mode_mid)
  html <- paste0(
    "<label style='", lbl_style, "'>", t$dil_mode, "</label>",
    "<select id='wrg_dil_mode' style='", sel_style, "margin-bottom:14px;' onchange='wrgUpdateDilemmas()'>",
    paste0(sprintf("<option value='%s'%s>%s</option>", modes, ifelse(modes == mode0, " selected", ""),
                   mode_lab[as.character(modes)]), collapse = ""), "</select>",
    "<label style='", lbl_style, "display:flex;justify-content:space-between;'><span>", t$dil_rmin,
    "</span><span id='wrg_dil_rmin_val' style='color:#8cc63f;'>", sprintf("%.2f", rmin0), "</span></label>",
    "<input id='wrg_dil_rmin' type='range' min='0' max='", length(rmins) - 1, "' step='1' value='",
    which(rmins == rmin0) - 1, "' style='width:100%;accent-color:#8cc63f;' oninput='wrgUpdateDilemmas()'>",
    "<div style='display:flex;justify-content:space-between;font-size:11px;color:#999;'><span>",
    sprintf("%.2f", min(rmins)), "</span><span>", sprintf("%.2f", max(rmins)), "</span></div>",
    "<label style='display:flex;align-items:center;gap:8px;margin-top:16px;font-size:13px;color:#555;cursor:pointer;'>",
    "<input id='wrg_dil_only' type='checkbox' style='accent-color:#8cc63f;' onchange='wrgUpdateDilemmas()'>",
    t$dil_toggle, "</label>")

  js <- paste0("
    var WRG_DIL = ", payload_json, ";
    var WRG_RMINS = ", jsonlite::toJSON(rmins), ";
    function wrgUpdateDilemmas() {
      var r = WRG_RMINS[parseInt(document.getElementById('wrg_dil_rmin').value, 10)];
      document.getElementById('wrg_dil_rmin_val').innerText = r.toFixed(2);
      var key = document.getElementById('wrg_dil_mode').value + '|' + r + '|' + (document.getElementById('wrg_dil_only').checked ? 'inv' : 'all');
      var v = WRG_DIL[key]; if (!v) return;
      var el = document.querySelector('#wrg_tab_1 .js-plotly-plot');
      if (el && window.Plotly) {
        Plotly.react(el, v.fig.data, v.fig.layout, {displayModeBar: false, responsive: true});
        el._wrgOrigAnnotations = null; // new data: forget the cached full-length labels
        wrgFitDilemmaPlot(el);
      }
      document.getElementById('wrg_dil_tbody').innerHTML = v.rows;
      ['n', 'pid', 'iid', 'picid'].forEach(function(k) { document.getElementById('wrg_kpi_' + k).innerText = v.kpi[k]; });
    }
  ")

  # The dilemma plot reserves large fixed pixel margins (see repgrid_dilemmas())
  # to fit the pole-label annotations either side of the plot. Those margins
  # don't shrink with the container, so on a narrow GridStack card they can
  # swallow the whole width and collapse the two marker columns onto each
  # other. wrgFitDilemmaPlot() re-scales the margins to the plot's actual
  # rendered width and trims label text with an ellipsis to whatever fits,
  # measured with canvas so it can never overflow the margin regardless of
  # font/margin choice. It re-runs on load, on resize (an iframe's own window
  # fires 'resize' when GridStack changes its size) and after every
  # Plotly.react() triggered by the settings menu. The full-length label text
  # is cached once per dataset (el._wrgOrigAnnotations) so enlarging the
  # widget again always re-expands from the original text instead of
  # re-trimming an already-shortened string.
  fit_js <- "
    function wrgMeasure(text, font) {
      var c = wrgMeasure._c || (wrgMeasure._c = document.createElement('canvas'));
      var ctx = c.getContext('2d');
      ctx.font = font;
      return ctx.measureText(text).width;
    }
    function wrgFitText(text, font, maxWidth) {
      if (maxWidth <= 0) return '';
      if (wrgMeasure(text, font) <= maxWidth) return text;
      var lo = 0, hi = text.length;
      while (lo < hi) {
        var mid = Math.ceil((lo + hi) / 2);
        var candidate = text.slice(0, mid) + '\\u2026';
        if (wrgMeasure(candidate, font) <= maxWidth) lo = mid; else hi = mid - 1;
      }
      return lo === 0 ? '' : text.slice(0, lo) + '\\u2026';
    }
    function wrgFitDilemmaPlot(el) {
      if (!el || !el.layout || !window.Plotly) return;
      // Base this on the widget's own document width, not el.clientWidth:
      // the shared responsive CSS already shrinks .html-widget by a fixed
      // gutter at narrow widths (room for the floating info/export
      // buttons), and computing margins from that ALREADY-shrunk width
      // double-counts the gutter, starving the actual plotting area far
      // below its minimum. Working from the full width and then forcing
      // el's own rendered width via JS (below) keeps the two independent.
      var totalW = document.body.clientWidth; if (!totalW) return;
      if (!el._wrgOrigAnnotations) el._wrgOrigAnnotations = JSON.parse(JSON.stringify(el.layout.annotations || []));
      var orig = el._wrgOrigAnnotations;

      // Mirror the two breakpoints .wt_responsive_css() uses to reserve
      // room for the floating info/export buttons, so we know up front how
      // much of totalW that gutter is already spoken for.
      var iconGutter = totalW <= 340 ? 38 : (totalW <= 520 ? 46 : 0);
      var avail = totalW - iconGutter;

      // Reserve a plausible fraction of the width for each side's labels,
      // but never let the inner plotting area (where the dots/lines live)
      // drop below a usable minimum - the content area's width always
      // wins over how much room the margins/labels get, and over the
      // icon-button gutter if it comes to that.
      var marginL = Math.min(310, totalW * 0.46);
      var marginR = Math.min(350, totalW * 0.50);
      var minInner = 140;
      if (avail - marginL - marginR < minInner) {
        var scale = Math.max(0, (avail - minInner) / (marginL + marginR));
        marginL *= scale; marginR *= scale;
      }
      marginL = Math.max(20, marginL); marginR = Math.max(20, marginR);
      if (avail - marginL - marginR < minInner) {
        // Margins are already at their floor and it's still not enough:
        // claim back whatever the icon gutter can spare.
        iconGutter = Math.max(0, totalW - marginL - minInner - marginR);
      }
      // El itself may be narrower than totalW (that same icon-button
      // gutter, applied via CSS); override it here so the plot always
      // renders at exactly the width this function just computed.
      var elWidth = totalW - iconGutter;
      el.style.setProperty('width', elWidth + 'px', 'important');
      var w = elWidth;
      var fontPx = Math.round(Math.max(9, Math.min(12, w / 60)));
      // Match Plotly's actual annotation font exactly (it falls back to
      // Verdana here since Open Sans isn't loaded in this iframe, and
      // Verdana is noticeably wider than a generic sans-serif) - otherwise
      // this measurement under-estimates the rendered width and the text
      // spills past its budget into the marker next to it.
      var font = 'bold ' + fontPx + 'px \"Open Sans\", verdana, arial, sans-serif';
      var xshiftPad = 20, safety = 26;

      var anns = orig.map(function(a) {
        var b = Object.assign({}, a);
        if (a.xanchor) {
          var m = /^<b>([\\s\\S]*?)<\\/b> - ([\\s\\S]*)$/.exec(a.text);
          var budget = (a.xanchor === 'right' ? marginL : marginR) - xshiftPad - safety;
          if (m) {
            var full = m[1] + ' - ' + m[2];
            var fitted = wrgFitText(full, font, budget);
            var sep = fitted.indexOf(' - ');
            b.text = sep >= 0 ? ('<b>' + fitted.slice(0, sep) + '</b>' + fitted.slice(sep)) : ('<b>' + fitted + '</b>');
          }
          b.font = Object.assign({}, a.font, {size: fontPx});
        } else if (a.bgcolor) {
          b.font = Object.assign({}, a.font, {size: Math.max(9, fontPx - 1)});
        }
        return b;
      });
      Plotly.relayout(el, {'margin.l': marginL, 'margin.r': marginR, annotations: anns});

      // Marker size follows the same scale as the label font so the dots
      // stay visually in proportion as the widget shrinks or grows; the
      // two marker traces (discrepant, congruent) are always the last two
      // in el.data, after one line trace per dilemma.
      var markerSize = Math.max(6, Math.min(11, Math.round(fontPx * 0.9167)));
      var n = el.data.length;
      if (n >= 2) Plotly.restyle(el, {'marker.size': markerSize}, [n - 2, n - 1]);
    }
    function wrgDilemmaPlotEl() { return document.querySelector('#wrg_tab_1 .js-plotly-plot'); }
    window.addEventListener('load', function() {
      setTimeout(function() { wrgFitDilemmaPlot(wrgDilemmaPlotEl()); }, 50);
      // Resizing the GridStack card resizes this iframe's CSS box, but that
      // does NOT fire a native 'resize' event inside the iframe's own
      // window (only an actual top-level window resize does), and a
      // ResizeObserver watching document.body from inside this same iframe
      // was unreliable in testing (it never notified for a resize driven
      // purely by the parent changing the iframe's box). Polling the
      // rendered width is cheap and has no such edge cases - it just
      // re-fits whenever the width actually changed since the last check.
      var wrgLastWidth = document.body.clientWidth;
      setInterval(function() {
        var w = document.body.clientWidth;
        if (w && w !== wrgLastWidth) {
          wrgLastWidth = w;
          wrgFitDilemmaPlot(wrgDilemmaPlotEl());
        }
      }, 300);
    });
    window.addEventListener('resize', function() {
      wrgFitDilemmaPlot(wrgDilemmaPlotEl());
    });
  "

  csv_js <- paste0(fit_js, sprintf("
    function downloadDilemmasCSV() {
      var rows = [['%s','%s','%s']];
      document.querySelectorAll('#wrg_dil_table tbody tr').forEach(function(r) {
        if (r.cells.length === 3) rows.push(Array.from(r.cells).map(function(c){ return '\"' + c.innerText.trim().replace(/\"/g, '\"\"') + '\"'; }));
      });
      var a = document.createElement('a');
      a.href = 'data:text/csv;charset=utf-8,' + encodeURIComponent(rows.map(function(r){ return r.join(','); }).join('\\n'));
      a.download = 'RepGrid_Implicative_Dilemmas.csv';
      a.click();
    }", t$congruent, t$discrepant, t$correlation))

  .rg_tabbed_widget(list(plot, htmltools::HTML(table_html)),
                    c(t$dilemma_graph_tab, t$dilemma_table_tab),
                    c("RepGrid_Dilemmas", "RepGrid_Dilemmas_Table"),
                    t$info_text_dilemma, t,
                    settings = list(html = html, js = js),
                    export_js = c(NA, "downloadDilemmasCSV()"),
                    extra_js = csv_js)
}


# widget_repgrid_indices ------------------------------------------------------

#' RepGrid Cognitive Indices Widget for Psychlab
#'
#' @description Creates a three-tab HTML widget: Tab 1 lists the global
#'   cognitive indices of the grid with a short description, Tab 2 and Tab 3
#'   give intensity, polarization and conflict for each construct and each
#'   element (see \code{\link{repgrid_indices}}).
#'
#' @param x An \code{OpenRepGrid} \code{repgrid} object. The self is assumed
#'   in the first column and the ideal in the last one.
#' @param lang Language for the UI. \code{"en"} (default) or \code{"es"}.
#'
#' @return A \code{browsable} HTML object with a three-tab interface.
#' @export
widget_repgrid_indices <- function(x, lang = "en") {
  if (!lang %in% c("en", "es")) lang <- "en"
  t <- wt_i18n(lang)
  ix <- repgrid_indices(x)
  esc <- function(v) htmltools::htmlEscape(v)
  fmt <- function(v, unit = "") {
    unit <- rep_len(unit, length(v))
    out <- ifelse(unit == "%", sprintf("%.2f%%", v), sprintf("%.2f", v))
    ifelse(is.na(v) | is.nan(v), "\u2014", out)
  }

  th <- function(label, tbl, i, right = FALSE) sprintf(
    "<th onclick=\"wrgSortTable('%s', %d)\" style='cursor:pointer;user-select:none;%s'>%s <span style='font-size:13px;color:#aaa;'>&#8597;</span></th>",
    tbl, i, if (right) "text-align:right;" else "", label)
  shell <- function(inner) paste0(
    "<div style='width:100%;height:100%;overflow:auto;padding:16px clamp(16px, 10vw, 70px) 60px 20px;box-sizing:border-box;font-family:Inter,Roboto,sans-serif;'>",
    inner, "</div>")

  # Tab 1: global indices, grouped
  g <- ix$global
  rows <- character(0)
  for (grp in unique(g$group)) {
    rows <- c(rows, sprintf("<tr class='wrg-grp'><td colspan='3'>%s</td></tr>", esc(t[[paste0("idx_g_", grp)]])))
    sub <- g[g$group == grp, ]
    rows <- c(rows, sprintf(
      "<tr><td style='font-weight:600;white-space:nowrap;'>%s</td><td style='text-align:right;font-weight:700;color:#333;'>%s</td><td style='color:#777;'>%s</td></tr>",
      esc(unlist(t[paste0("idx_n_", sub$key)])), fmt(sub$value, sub$unit), esc(unlist(t[paste0("idx_d_", sub$key)]))))
  }
  tab1 <- shell(paste0(
    "<table id='wrg_idx_global' class='wrg-idx'><thead><tr><th>", t$idx_index,
    "</th><th style='text-align:right;'>", t$idx_value, "</th><th>", t$idx_desc, "</th></tr></thead><tbody>",
    paste(rows, collapse = "\n"), "</tbody></table>"))

  # Tabs 2 and 3: per construct / per element, sortable
  per_table <- function(df, id, first_label, first_col) {
    body <- paste0(sprintf(
      "<tr><td>%s</td><td style='text-align:right;'>%s</td><td style='text-align:right;'>%s</td><td style='text-align:right;'>%s</td></tr>",
      esc(df[[first_col]]), fmt(df$intensity), fmt(df$polarization, "%"), fmt(df$conflict, "%")), collapse = "\n")
    shell(paste0(
      "<table id='", id, "' class='wrg-idx'><thead><tr>", th(first_label, id, 0),
      th(t$idx_intensity, id, 1, TRUE), th(t$idx_polarization, id, 2, TRUE), th(t$idx_conflict, id, 3, TRUE),
      "</tr></thead><tbody>", body, "</tbody></table>"))
  }
  tab2 <- per_table(ix$constructs, "wrg_idx_constructs", t$idx_construct, "construct")
  tab3 <- per_table(ix$elements, "wrg_idx_elements", t$idx_element, "element")

  css <- "<style>
    .wrg-idx{width:100%;border-collapse:collapse;font-size:13px;}
    .wrg-idx thead th{padding:10px 14px;text-align:left;font-weight:600;color:#555;font-size:12px;background:#f8f9fa;border-bottom:2px solid #8cc63f;white-space:nowrap;}
    .wrg-idx thead th:hover{background:#eef7e0;}
    .wrg-idx tbody tr{border-bottom:1px solid #f0f0f0;}
    .wrg-idx tbody tr:hover{background:#f4f9ef;}
    .wrg-idx tbody td{padding:9px 14px;color:#333;}
    .wrg-idx tr.wrg-grp td{background:#fbfdf9;color:#8cc63f;font-weight:700;font-size:12px;text-transform:uppercase;letter-spacing:.4px;padding:12px 14px 6px;border-bottom:1px solid #e3efd0;}
  </style>"

  js <- "
    var _wrgDir = {};
    function wrgSortTable(id, col) {
      var tb = document.querySelector('#' + id + ' tbody');
      var rows = Array.from(tb.querySelectorAll('tr'));
      var k = id + col; _wrgDir[k] = !_wrgDir[k]; var asc = _wrgDir[k];
      rows.sort(function(a, b) {
        var va = a.cells[col].innerText.trim(), vb = b.cells[col].innerText.trim();
        var na = parseFloat(va), nb = parseFloat(vb);
        if (!isNaN(na) && !isNaN(nb)) return asc ? na - nb : nb - na;
        return asc ? va.localeCompare(vb) : vb.localeCompare(va);
      });
      rows.forEach(function(r) { tb.appendChild(r); });
    }
    function downloadIndicesCSV(id, name) {
      var rows = [];
      document.querySelectorAll('#' + id + ' tr').forEach(function(r) {
        if (r.style.display === 'none') return;
        rows.push(Array.from(r.cells).map(function(c) { return '\"' + c.innerText.replace(/[\\u2195\\u2191\\u2193]/g, '').trim().replace(/\"/g, '\"\"') + '\"'; }).join(','));
      });
      var a = document.createElement('a');
      a.href = 'data:text/csv;charset=utf-8,' + encodeURIComponent(rows.join('\\n'));
      a.download = name + '.csv';
      a.click();
    }"

  .rg_tabbed_widget(
    list(htmltools::HTML(paste0(css, tab1)), htmltools::HTML(tab2), htmltools::HTML(tab3)),
    c(t$idx_tab_global, t$idx_tab_constructs, t$idx_tab_elements),
    c("RepGrid_Indices", "RepGrid_Indices_Constructs", "RepGrid_Indices_Elements"),
    t$info_text_indices, t,
    export_js = c("downloadIndicesCSV('wrg_idx_global', 'RepGrid_Indices')",
                  "downloadIndicesCSV('wrg_idx_constructs', 'RepGrid_Indices_Constructs')",
                  "downloadIndicesCSV('wrg_idx_elements', 'RepGrid_Indices_Elements')"),
    extra_js = js)
}
