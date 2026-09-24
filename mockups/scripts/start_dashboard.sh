#!/bin/bash
# Regenerates the R-rendered widgets/*.html files and serves the mockups/ folder.
# Self-locating: works no matter where it's invoked from.
cd "$(dirname "$0")/.." || exit 1

echo "Generando widgets desde R..."
Rscript --vanilla -e "
devtools::load_all('..')
library(htmltools)
data(su_wimp)
data(example_wimp)
save_html(widget_digraph(su_wimp, lang='es'), 'widgets/widget_digraph.html')
save_html(widget_centrality(su_wimp, lang='es'), 'widgets/widget_centrality.html')
save_html(widget_implications(su_wimp, lang='es'), 'widgets/widget_implications.html')
save_html(widget_simulation(su_wimp, lang='es'), 'widgets/widget_simulation.html')
save_html(widget_adjustment(example_wimp, lang='es'), 'widgets/widget_adjustment_single.html')
save_html(widget_adjustment(example_wimp, example_wimp_post, lang='es'), 'widgets/widget_adjustment_compare.html')
save_html(widget_repgrid_biplot(OpenRepGrid::feixas2004, lang='es'), 'widgets/widget_repgrid_biplot.html')
save_html(widget_repgrid_cluster(OpenRepGrid::feixas2004, lang='es'), 'widgets/widget_repgrid_cluster.html')
save_html(widget_repgrid_dilemmas(OpenRepGrid::feixas2004, lang='es'), 'widgets/widget_repgrid_dilemmas.html')
save_html(widget_repgrid_indices(OpenRepGrid::feixas2004, lang='es'), 'widgets/widget_repgrid_indices.html')
"

if [ $? -eq 0 ]; then
    echo "¡Widgets generados! Iniciando servidor local en http://localhost:8080/dashboard_test.html"
    # Try to open browser (Linux)
    xdg-open "http://localhost:8080/dashboard_test.html" 2>/dev/null &
    python3 -m http.server 8080
else
    echo "Error generando los widgets."
fi
