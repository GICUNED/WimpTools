# WimpTools Localization (i18n)
# ============================================================
# Internal helper function. NOT exported.
# Usage: t <- wt_i18n("es")   # returns the Spanish translation list
# To add a new language, add a new entry to the `dict` list below.
# ============================================================

wt_i18n <- function(lang = "en") {

  dict <- list(

    # ── English ────────────────────────────────────────────────
    en = list(
      # Tab labels
      self_digraph   = "Self Digraph",
      weight_matrix  = "Weight Matrix",

      # Floating buttons (tooltips)
      export_png     = "Export PNG",
      info           = "Info",
      fullscreen     = "Fullscreen",

      # Area / grouping panel
      areas          = "Areas",
      group_by       = "Group by",
      none           = "None",
      cluster_areas  = "Cluster Areas",

      # Visualization panel — labels
      vis_options        = "Settings",
      color_palette      = "Color Palette",
      layout_algo        = "Layout Algorithm",
      edge_filter        = "Edge Filter",
      edge_opacity       = "Edge Opacity",
      hide_direct        = "Hide Direct (Positive)",
      node_spacing       = "Node Spacing",
      node_size          = "Node Size",
      text_size          = "Text Size",
      eraser_mode        = "Eraser Mode",
      reset              = "Reset",
      visible_constructs = "Visible Constructs",
      all                = "All",
      none_btn           = "None",
      congruence         = "Congruence",
      category           = "Category",

      # Simulation panel — labels
      sim_settings   = "Settings",
      thr_function   = "Threshold Function",
      saturation     = "Saturation",
      tanh           = "Hyperbolic (Tanh)",
      linear         = "Linear (No limit)",
      sim_depth      = "Simulation Depth (Iter)",
      playback_speed = "Playback Speed",
      scenario       = "Scenario",
      reset_scenario = "Reset Scenario",
      network_view   = "Network View",
      pcsd_chart     = "PCSD Chart",

      # Heatmap settings
      hm_settings    = "Settings",
      sort_by        = "Sort By",
      sort_original  = "Original",
      sort_weight    = "Absolute Weight",
      sort_connect   = "Connectivity",
      show_density   = "Show Density \u03c1(G)",
      show_values    = "Show Cell Values",
      filter_constructs = "Filter Constructs",
      pal_redgreen   = "Red - Green",
      pal_redblue    = "Red - Blue",
      pal_orangepurple = "Orange - Purple",
      pal_greyscale  = "Greyscale",

      # Info modal body text
      info_text_digraph  = paste0(
        "The Self Digraph represents the structure of the cognitive system. ",
        "Nodes correspond to personal constructs, and edges (arrows) indicate ",
        "directions of implication or influence. Larger or central nodes have a ",
        "greater impact on the system, and their color reflects their alignment ",
        "with the ideal pole."
      ),
      info_text_heatmap  = paste0(
        "The Weight Matrix shows the quantitative strength of influences between ",
        "constructs. Rows represent the influence capacity of a construct (source) ",
        "on others, while columns represent the degree to which a construct is ",
        "affected (destination). More intense colors indicate a stronger structural ",
        "relationship."
      ),
      info_text_ideal_digraph = paste0(
        "The Ideal Digraph visualizes the cognitive system as it would be if all constructs ",
        "were perfectly aligned with the individual's ideal pole. Nodes represent constructs, ",
        "and edges represent implications. Exploring this structure reveals the desired or ",
        "aspirational state of the psychological system."
      ),
      info_text_if = paste0(
        "The Impact and Feedback (IF) Plot evaluates each construct based on two dimensions: ",
        "Impact (ability to influence other constructs) and Feedback (susceptibility to be influenced). ",
        "This helps identify which constructs drive the system and which ones are more reactive."
      ),
      info_text_hypo = paste0(
        "The Hypothetical Scenarios visualization allows you to explore what would happen ",
        "if specific constructs were modified. It simulates the propagation of change through ",
        "the implication network, showing the systemic consequences of local interventions."
      ),

      # Dropdown option labels — color palettes
      palette_labels = list(
        `red/green`  = "Red / Green",
        `grey scale` = "Grey Scale",
        colorblind   = "Colorblind",
        pastel       = "Pastel",
        dark         = "Dark",
        viridis      = "Viridis"
      ),

      # Dropdown option labels — layouts
      layout_labels = list(
        graphopt = "Graphopt",
        circle   = "Circle",
        mds      = "MDS",
        grid     = "Grid",
        tree     = "Tree",
        areas    = "Areas"
      ),

      # Centrality widget — tab labels
      pb_plot_tab    = "PB Plot",
      centrality_tab = "Centrality Indices",

      # Centrality widget — table column headers
      construct     = "Construct",
      presence      = "Presence (P)",
      balance       = "Balance (B)",
      degree        = "Degree",
      closeness     = "Closeness",
      betweenness   = "Betweenness",

      # Centrality widget — PB plot axes
      pb_axis_x     = "PRESENCE",
      pb_axis_y     = "BALANCE",

      # Centrality widget — info modal texts
      info_text_pb = paste0(
        "The PB Plot represents each construct in a two-dimensional space defined ",
        "by Presence (P) and Balance (B). Presence reflects how connected a ",
        "construct is to the rest of the system (high P = high influence). ",
        "Balance indicates whether a construct predominantly influences others ",
        "(positive B) or is influenced by them (negative B)."
      ),
      info_text_centrality = paste0(
        "This table displays key network metrics for each construct. ",
        "Degree measures the total weight of connections. Closeness reflects ",
        "how quickly a construct can reach the rest of the system. ",
        "Betweenness identifies constructs that act as bridges in the ",
        "shortest paths between other constructs."
      ),
      info_text_sim_network = paste0(
        "The Simulated Digraph visually represents the evolving state of the cognitive system. ",
        "Use the left panel to define a hypothetical scenario by altering the initial activation ",
        "of specific constructs. The network will then recursively update, allowing you to observe ",
        "how systemic dependencies propagate changes and stabilize over time."
      ),
      info_text_sim_pcsd = paste0(
        "The Personal Construct System Dynamics (PCSD) Chart tracks the temporal trajectory of each ",
        "construct throughout the simulation. The Y-axis represents the self-differential state, ",
        "ranging from -2 to 2, illustrating the system's structural constraints and conflict resolution paths."
      )
    ),

    # ── Español ────────────────────────────────────────────────
    es = list(
      # Etiquetas de pestañas
      self_digraph   = "Digrafo del Self",
      weight_matrix  = "Matriz de Pesos",

      # Botones flotantes (tooltips)
      export_png     = "Exportar PNG",
      info           = "Info",
      fullscreen     = "Pantalla Completa",

      # Panel de áreas / agrupación
      areas          = "\u00c1reas",
      group_by       = "Agrupar por",
      none           = "Ninguno",
      cluster_areas  = "Agrupar \u00c1reas",

      # Panel de visualización — etiquetas
      vis_options        = "Ajustes",
      color_palette      = "Paleta de Color",
      layout_algo        = "Algoritmo de Dise\u00f1o",
      edge_filter        = "Filtro de Aristas",
      edge_opacity       = "Opacidad Aristas",
      hide_direct        = "Ocultar Directas (Positivas)",
      node_spacing       = "Espaciado de Nodos",
      node_size          = "Tama\u00f1o de Nodo",
      text_size          = "Tama\u00f1o de Texto",
      eraser_mode        = "Modo Borrador",
      reset              = "Restaurar",
      visible_constructs = "Constructos Visibles",
      all                = "Todos",
      none_btn           = "Ninguno",
      congruence         = "Congruencia",
      category           = "Categor\u00eda",

      # Panel de simulación — etiquetas
      sim_settings   = "Ajustes",
      thr_function   = "Funci\u00f3n de Activaci\u00f3n",
      saturation     = "Saturaci\u00f3n",
      tanh           = "Hiperb\u00f3lica (Tanh)",
      linear         = "Lineal (Sin l\u00edmite)",
      sim_depth      = "Profundidad (Iter)",
      playback_speed = "Velocidad de Repr.",
      scenario       = "Escenario",
      reset_scenario = "Reiniciar Escenario",
      network_view   = "Vista de Red",
      pcsd_chart     = "Gr\u00e1fico PCSD",

      # Heatmap settings
      hm_settings    = "Ajustes",
      sort_by        = "Ordenar por",
      sort_original  = "Original",
      sort_weight    = "Peso Absoluto",
      sort_connect   = "Conectividad",
      show_density   = "Mostrar Densidad \u03c1(G)",
      show_values    = "Mostrar Valores",
      filter_constructs = "Filtrar Constructos",
      pal_redgreen   = "Rojo - Verde",
      pal_redblue    = "Rojo - Azul",
      pal_orangepurple = "Naranja - Morado",
      pal_greyscale  = "Escala de grises",

      # Texto de los modales de información
      info_text_digraph  = paste0(
        "El Digrafo del Self representa la estructura del sistema cognitivo. ",
        "Los nodos corresponden a los constructos personales, y las aristas (flechas) ",
        "indican las direcciones de implicaci\u00f3n o influencia. Los nodos m\u00e1s grandes o ",
        "centrales tienen mayor impacto en el sistema, y su color refleja su alineaci\u00f3n ",
        "con el polo ideal."
      ),
      info_text_heatmap  = paste0(
        "La Matriz de Pesos (Weight Matrix) muestra la fuerza cuantitativa de las ",
        "influencias entre los constructos. Las filas representan la capacidad de ",
        "influencia de un constructo (origen) sobre los dem\u00e1s, mientras que las columnas ",
        "representan el grado en que un constructo es afectado (destino). Colores m\u00e1s ",
        "intensos indican una relaci\u00f3n estructural m\u00e1s fuerte."
      ),
      info_text_ideal_digraph = paste0(
        "El Digrafo del Ideal visualiza el sistema cognitivo tal como sería si todos los constructos ",
        "estuvieran perfectamente alineados con el polo ideal del individuo. Los nodos representan ",
        "constructos y las aristas implicaciones. Explorar esta estructura revela el estado ",
        "deseado o aspiracional del sistema psicológico."
      ),
      info_text_if = paste0(
        "El Gráfico de Impacto y Feedback evalúa cada constructo en base a dos dimensiones: ",
        "Impacto (capacidad de influir sobre otros constructos) y Feedback (susceptibilidad de ser influido). ",
        "Esto permite identificar qué constructos dirigen el sistema y cuáles son más reactivos."
      ),
      info_text_hypo = paste0(
        "La visualización de Escenarios Hipotéticos permite explorar qué sucedería si se ",
        "modificaran constructos específicos. Simula la propagación del cambio a través de la red ",
        "de implicaciones, mostrando las consecuencias sistémicas de intervenciones locales."
      ),

      # Etiquetas de opciones desplegables — paletas de color
      palette_labels = list(
        `red/green`  = "Rojo / Verde",
        `grey scale` = "Escala de Grises",
        colorblind   = "Daltonismo",
        pastel       = "Pastel",
        dark         = "Oscuro",
        viridis      = "Viridis"
      ),

      # Etiquetas de opciones desplegables — algoritmos de disposición
      layout_labels = list(
        graphopt = "Graphopt",
        circle   = "Circular",
        mds      = "MDS",
        grid     = "Cuadr\u00edcula",
        tree     = "\u00c1rbol",
        areas    = "\u00c1reas"
      ),

      # Widget de centralidad — etiquetas de pestañas
      pb_plot_tab    = "Gr\u00e1fico PB",
      centrality_tab = "\u00cdndices de centralidad",

      # Widget de centralidad — cabeceras de la tabla
      construct     = "Constructo",
      presence      = "Presencia (P)",
      balance       = "Balance (B)",
      degree        = "Grado",
      closeness     = "Cercan\u00eda",
      betweenness   = "Intermediaci\u00f3n",

      # Widget de centralidad — ejes del gráfico PB
      pb_axis_x     = "PRESENCIA",
      pb_axis_y     = "BALANCE",

      # Widget de centralidad — textos de los modales de información
      info_text_pb = paste0(
        "El Gr\u00e1fico PB representa cada constructo en un espacio bidimensional ",
        "definido por la Presencia (P) y el Balance (B). La Presencia refleja ",
        "c\u00f3mo de conectado est\u00e1 un constructo con el resto del sistema ",
        "(P alto = alta influencia). El Balance indica si un constructo ",
        "predominantemente influye sobre los dem\u00e1s (B positivo) o es influido ",
        "por ellos (B negativo)."
      ),
      info_text_centrality = paste0(
        "Esta tabla muestra los principales \u00edndices de centralidad para cada ",
        "constructo. El Grado mide la conectividad total (suma de pesos). ",
        "La Cercan\u00eda refleja con qu\u00e9 rapidez un constructo puede alcanzar al ",
        "resto del sistema. La Intermediaci\u00f3n cuenta cu\u00e1ntas veces un constructo ",
        "act\u00faa como puente en los caminos m\u00e1s cortos entre otros constructos."
      ),
      info_text_sim_network = paste0(
        "El Digrafo Simulado representa visualmente la evoluci\u00f3n del estado del sistema cognitivo. ",
        "Utiliza el panel izquierdo para definir un escenario hipot\u00e9tico alterando la activaci\u00f3n inicial ",
        "de constructos espec\u00edficos. La red se actualizar\u00e1 recursivamente, permitiendo observar c\u00f3mo ",
        "las dependencias sist\u00e9micas propagan los cambios y se estabilizan con el tiempo."
      ),
      info_text_sim_pcsd = paste0(
        "El Gr\u00e1fico de Din\u00e1mica de Sistemas de Constructos Personales (PCSD) rastrea la ",
        "trayectoria temporal de cada constructo a lo largo de la simulaci\u00f3n. El eje Y representa ",
        "el estado diferencial, variando de -2 a 2, ilustrando las restricciones estructurales ",
        "y las v\u00edas de resoluci\u00f3n de conflictos del sistema."
      )
    )
  )

  if (!lang %in% names(dict)) lang <- "en"
  dict[[lang]]
}
