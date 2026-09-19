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
      biplot_2d_tab  = "Biplot 2D",
      biplot_3d_tab  = "Biplot 3D",
      cluster_constructs_tab = "Constructs Cluster",
      cluster_elements_tab   = "Elements Cluster",
      dilemma_graph_tab = "Dilemmas Graph",
      dilemma_table_tab = "Dilemmas Table",
      dil_mode       = "Congruence criterion",
      dil_mode_diff  = "Self-ideal difference",
      dil_mode_mid   = "Polarity",
      dil_rmin       = "Minimum correlation",
      dil_toggle     = "Only show constructs with dilemmas",
      congruent      = "Congruent construct",
      discrepant     = "Discrepant construct",
      correlation    = "Correlation",
      n_dilemmas     = "Implicative dilemmas",
      no_dilemmas    = "No implicative dilemmas found",
      export_csv     = "Export CSV",
      info_text_dilemma = paste0(
        "An implicative dilemma appears when a construct that is congruent (the person is already where they want to be) ",
        "is positively correlated with a construct that is discrepant (the person wants to change), so changing one ",
        "implies changing the other. Lines join the pairs; thicker lines mean stronger correlation. ",
        "PID: percentage of dilemmas; IID: intensity; PICID: proportion of the intensity of the constructs involved."
      ),
      idx_tab_global = "Indices",
      idx_tab_constructs = "Constructs",
      idx_tab_elements = "Elements",
      idx_index = "Index", idx_value = "Value", idx_desc = "Description",
      idx_construct = "Construct", idx_element = "Element",
      idx_intensity = "Intensity", idx_polarization = "Polarization", idx_conflict = "Conflict",
      info_text_indices = paste0(
        "Cognitive indices of the repertory grid computed with OpenRepGrid. The first tab summarises the global indices ",
        "(structure, conflict and self-related); the other tabs give intensity, polarization and conflict for each construct ",
        "and each element. Click a column header to sort. The self is taken from the first column and the ideal from the last."
      ),
      idx_g_structure = "Cognitive structure", idx_g_conflict = "Conflict", idx_g_self = "Self and ideal",
      idx_n_int_total = "Intensity (total)", idx_d_int_total = "Overall degree of linkage among constructs and elements (Bannister). Higher values mean a more tightly organised system.",
      idx_n_int_c = "Intensity (constructs)", idx_d_int_c = "Mean intensity of the constructs.",
      idx_n_int_e = "Intensity (elements)", idx_d_int_e = "Mean intensity of the elements.",
      idx_n_pvaff = "PVAFF", idx_d_pvaff = "Percentage of variance accounted for by the first factor. A measure of cognitive complexity: high values indicate a simple, one-dimensional system.",
      idx_n_bieri = "Bieri", idx_d_bieri = "Proportion of matching ratings between pairs of constructs. Higher values mean lower cognitive differentiation.",
      idx_n_variability = "Variability", idx_d_variability = "Tendency of the ratings to gravitate towards both ends of the scale (Slater).",
      idx_n_bias = "Bias", idx_d_bias = "Tendency of the ratings to accumulate at one end of the scale (Slater).",
      idx_n_polarization = "Polarization", idx_d_polarization = "Percentage of ratings placed at the extreme points of the scale.",
      idx_n_conf1 = "Conflict (Slade & Sheehan)", idx_d_conf1 = "Percentage of imbalanced triads of constructs (balance theory).",
      idx_n_conf2 = "Conflict (Bassler et al.)", idx_d_conf2 = "Percentage of imbalanced triads, also taking into account the magnitude of the correlations.",
      idx_n_conf3 = "Conflict (Bell)", idx_d_conf3 = "Percentage of inconsistent element-construct-construct triads, based on distances.",
      idx_n_self_ideal = "Self - Ideal", idx_d_self_ideal = "Normalised distance between the self and the ideal self. Lower means closer.",
      idx_n_self_others = "Self - Others", idx_d_self_others = "Normalised distance between the self and the other elements.",
      idx_n_ideal_others = "Ideal - Others", idx_d_ideal_others = "Normalised distance between the ideal self and the other elements.",
      idx_n_dilemmatic = "Dilemmatic constructs", idx_d_dilemmatic = "Percentage of constructs where the ideal is rated at the mid-point of the scale.",
      cluster_dist   = "Distance measure",
      cluster_method = "Linkage method",
      info_text_cluster = paste0(
        "Hierarchical cluster analysis of the repertory grid. The dendrogram joins the most similar constructs (or elements) ",
        "first; the shorter the horizontal distance at which two branches merge, the more similar they are. ",
        "Hover a leaf to see the full label."
      ),
      info_text_biplot = paste0(
        "The biplot is a principal component representation of the repertory grid. Elements are shown as points ",
        "and constructs as vectors: elements lying in the direction of a construct vector are rated high on that pole. ",
        "The 3D tab adds the third component; drag to rotate."
      ),

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
      biplot_2d_tab  = "Biplot 2D",
      biplot_3d_tab  = "Biplot 3D",
      cluster_constructs_tab = "Cl\u00faster de Constructos",
      cluster_elements_tab   = "Cl\u00faster de Elementos",
      dilemma_graph_tab = "Gr\u00e1fico de Dilemas",
      dilemma_table_tab = "Tabla de Dilemas",
      dil_mode       = "Criterio de congruencia",
      dil_mode_diff  = "Diferencia self-ideal",
      dil_mode_mid   = "Polaridad",
      dil_rmin       = "Correlaci\u00f3n m\u00ednima",
      dil_toggle     = "Mostrar solo constructos con dilemas",
      congruent      = "Constructo congruente",
      discrepant     = "Constructo discrepante",
      correlation    = "Correlaci\u00f3n",
      n_dilemmas     = "Dilemas implicativos",
      no_dilemmas    = "No se han encontrado dilemas implicativos",
      export_csv     = "Exportar CSV",
      info_text_dilemma = paste0(
        "Un dilema implicativo aparece cuando un constructo congruente (la persona ya est\u00e1 donde quiere estar) ",
        "correlaciona positivamente con un constructo discrepante (la persona quiere cambiar), de modo que cambiar uno ",
        "implica cambiar el otro. Las l\u00edneas unen los pares; cuanto m\u00e1s gruesas, mayor correlaci\u00f3n. ",
        "PID: porcentaje de dilemas; IID: intensidad; PICID: proporci\u00f3n de la intensidad de los constructos implicados."
      ),
      idx_tab_global = "\u00cdndices",
      idx_tab_constructs = "Constructos",
      idx_tab_elements = "Elementos",
      idx_index = "\u00cdndice", idx_value = "Valor", idx_desc = "Descripci\u00f3n",
      idx_construct = "Constructo", idx_element = "Elemento",
      idx_intensity = "Intensidad", idx_polarization = "Polarizaci\u00f3n", idx_conflict = "Conflicto",
      info_text_indices = paste0(
        "\u00cdndices cognitivos de la rejilla de repertorio calculados con OpenRepGrid. La primera pesta\u00f1a resume los ",
        "\u00edndices globales (estructura, conflicto y relativos al yo); las otras dan la intensidad, la polarizaci\u00f3n y el ",
        "conflicto de cada constructo y de cada elemento. Pulsa una cabecera para ordenar. El yo se toma de la primera columna y el ideal de la \u00faltima."
      ),
      idx_g_structure = "Estructura cognitiva", idx_g_conflict = "Conflicto", idx_g_self = "Yo e ideal",
      idx_n_int_total = "Intensidad (total)", idx_d_int_total = "Grado global de vinculaci\u00f3n entre constructos y elementos (Bannister). Valores altos indican un sistema m\u00e1s organizado.",
      idx_n_int_c = "Intensidad (constructos)", idx_d_int_c = "Intensidad media de los constructos.",
      idx_n_int_e = "Intensidad (elementos)", idx_d_int_e = "Intensidad media de los elementos.",
      idx_n_pvaff = "PVAFF", idx_d_pvaff = "Porcentaje de varianza explicada por el primer factor. Medida de complejidad cognitiva: valores altos indican un sistema simple, unidimensional.",
      idx_n_bieri = "Bieri", idx_d_bieri = "Proporci\u00f3n de coincidencias en las puntuaciones entre pares de constructos. Valores altos indican menor diferenciaci\u00f3n cognitiva.",
      idx_n_variability = "Variabilidad", idx_d_variability = "Tendencia de las puntuaciones a situarse hacia ambos extremos de la escala (Slater).",
      idx_n_bias = "Sesgo", idx_d_bias = "Tendencia de las puntuaciones a acumularse en un extremo de la escala (Slater).",
      idx_n_polarization = "Polarizaci\u00f3n", idx_d_polarization = "Porcentaje de puntuaciones situadas en los puntos extremos de la escala.",
      idx_n_conf1 = "Conflicto (Slade y Sheehan)", idx_d_conf1 = "Porcentaje de tr\u00edadas de constructos desequilibradas (teor\u00eda del equilibrio).",
      idx_n_conf2 = "Conflicto (Bassler et al.)", idx_d_conf2 = "Porcentaje de tr\u00edadas desequilibradas, teniendo en cuenta tambi\u00e9n la magnitud de las correlaciones.",
      idx_n_conf3 = "Conflicto (Bell)", idx_d_conf3 = "Porcentaje de tr\u00edadas elemento-constructo-constructo inconsistentes, basado en distancias.",
      idx_n_self_ideal = "Yo - Ideal", idx_d_self_ideal = "Distancia normalizada entre el yo y el yo ideal. Menor valor, mayor cercan\u00eda.",
      idx_n_self_others = "Yo - Otros", idx_d_self_others = "Distancia normalizada entre el yo y el resto de elementos.",
      idx_n_ideal_others = "Ideal - Otros", idx_d_ideal_others = "Distancia normalizada entre el yo ideal y el resto de elementos.",
      idx_n_dilemmatic = "Constructos dilem\u00e1ticos", idx_d_dilemmatic = "Porcentaje de constructos en los que el ideal se punt\u00faa en el punto medio de la escala.",
      cluster_dist   = "Medida de distancia",
      cluster_method = "M\u00e9todo de enlace",
      info_text_cluster = paste0(
        "An\u00e1lisis de cl\u00faster jer\u00e1rquico de la rejilla de repertorio. El dendrograma une primero los constructos ",
        "(o elementos) m\u00e1s parecidos; cuanto menor es la distancia horizontal a la que se unen dos ramas, m\u00e1s similares son. ",
        "Pasa el rat\u00f3n por una hoja para ver la etiqueta completa."
      ),
      info_text_biplot = paste0(
        "El biplot es una representaci\u00f3n en componentes principales de la rejilla de repertorio. Los elementos se ",
        "muestran como puntos y los constructos como vectores: los elementos situados en la direcci\u00f3n de un vector ",
        "puntúan alto en ese polo. La pesta\u00f1a 3D a\u00f1ade el tercer componente; arrastra para rotar."
      ),

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
