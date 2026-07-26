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
            <div style='background-color: rgba(255, 255, 255, 0.95); width: 32px; height: 32px; border-radius: 6px; box-shadow: 0 2px 10px rgba(0,0,0,0.1); border: 1px solid #ddd; display: flex; align-items: center; justify-content: center; cursor: pointer; transition: all 0.2s;' title='", t$info, "' onmouseover=\"this.style.backgroundColor='#f5f5f5'\" onmouseout=\"this.style.backgroundColor='rgba(255, 255, 255, 0.95)'\" onclick=\"var m=document.getElementById('heatmap_info_modal'); m.style.display=(m.style.display==='block'?'none':'block');\">
              <svg width='18' height='18' viewBox='0 0 24 24' fill='none' stroke='#333' stroke-width='2.5' stroke-linecap='round' stroke-linejoin='round'><circle cx='12' cy='12' r='10'></circle><line x1='12' y1='16' x2='12' y2='12'></line><line x1='12' y1='8' x2='12.01' y2='8'></line></svg>
            </div>
            <div id='hm_settings_placeholder'></div>
            <div style='background-color: rgba(255, 255, 255, 0.95); width: 32px; height: 32px; border-radius: 6px; box-shadow: 0 2px 10px rgba(0,0,0,0.1); border: 1px solid #ddd; display: flex; align-items: center; justify-content: center; cursor: pointer; transition: all 0.2s;' title='", t$export_png, "' onmouseover=\"this.style.backgroundColor='#f5f5f5'\" onmouseout=\"this.style.backgroundColor='rgba(255, 255, 255, 0.95)'\" onclick=\"downloadHeatmap()\">
              <svg width='18' height='18' viewBox='0 0 24 24' fill='none' stroke='#333' stroke-width='2.5' stroke-linecap='round' stroke-linejoin='round'><path d='M21 15v4a2 2 0 0 1-2 2H5a2 2 0 0 1-2-2v-4'></path><polyline points='7 10 12 15 17 10'></polyline><line x1='12' y1='15' x2='12' y2='3'></line></svg>
            </div>
            <div style='background-color: rgba(255, 255, 255, 0.95); width: 32px; height: 32px; border-radius: 6px; box-shadow: 0 2px 10px rgba(0,0,0,0.1); border: 1px solid #ddd; display: flex; align-items: center; justify-content: center; cursor: pointer; transition: all 0.2s;' title='", t$fullscreen, "' onmouseover=\"this.style.backgroundColor='#f5f5f5'\" onmouseout=\"this.style.backgroundColor='rgba(255, 255, 255, 0.95)'\" onclick=\"var el=this.closest('.wt-tab-container')||this.closest('.wt-tab-content'); if(!document.fullscreenElement){el.requestFullscreen().catch(e=>console.log(e))}else{document.exitFullscreen()}\">
              <svg width='18' height='18' viewBox='0 0 24 24' fill='none' stroke='#333' stroke-width='2.5' stroke-linecap='round' stroke-linejoin='round'><path d='M8 3H5a2 2 0 0 0-2 2v3m18 0V5a2 2 0 0 0-2-2h-3m0 18h3a2 2 0 0 0 2-2v-3M3 16v3a2 2 0 0 0 2 2h3'></path></svg>
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
    "<div style='flex:1; overflow:auto; padding:0 70px 60px 20px;'>",
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
      width: 320px;
      min-width: 320px;
      height: 100%;
      background: #f8f9fa;
      border-right: 1px solid #ddd;
      padding: 15px;
      overflow-y: hidden;
      overflow-x: hidden;
      display: flex;
      flex-direction: column;
    }
    .wsim-main {
      flex: 1;
      display: flex;
      flex-direction: column;
      height: 100%;
      overflow: hidden;
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
  plot_if <- if_barchart(x, ...)
  plot_if$width <- "100%"
  plot_if$height <- "100%"
  plot_if$sizingPolicy$defaultWidth <- "100%"
  plot_if$sizingPolicy$defaultHeight <- "100%"
  
  # Hypo Plot
  plot_hypo <- hypo_plot(x, ...)
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
