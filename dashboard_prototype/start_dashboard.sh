#!/bin/bash
echo "Generando widgets desde R..."
Rscript --vanilla -e "
devtools::load_all('..')
library(htmltools)
data(su_wimp)
save_html(widget_digraph(su_wimp, lang='es'), 'widget_digraph.html')
save_html(widget_centrality(su_wimp, lang='es'), 'widget_centrality.html')
save_html(widget_implications(su_wimp, lang='es'), 'widget_implications.html')
save_html(widget_simulation(su_wimp, lang='es'), 'widget_simulation.html')
"

if [ $? -eq 0 ]; then
    echo "¡Widgets generados! Iniciando servidor local en http://localhost:8080/dashboard_test.html"
    # Try to open browser (Linux)
    xdg-open "http://localhost:8080/dashboard_test.html" 2>/dev/null &
    python3 -m http.server 8080
else
    echo "Error generando los widgets."
fi
