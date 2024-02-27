## EMOTIONAL ADJUSTMENT ##

Adjustment_self_ideal <- function(wimp) {
  # Extraer los vectores ideal y self de la variable wimp
  vector_ideal <- wimp$ideal[[2]]
  vector_self <- wimp$self[[2]]

  # Calcular la distancia entre los vectores
  distancia <- sqrt(sum((vector_ideal - vector_self)^2))

  # Normalizar la distancia
  distancia_normalizada <- distancia / (2 * sqrt(length(vector_ideal)))

  # Calcular el coseno del ángulo entre los vectores (alineación)
  alineacion <- sum(vector_ideal * vector_self) / (sqrt(sum(vector_ideal^2)) * sqrt(sum(vector_self^2)))

  # Crear el vector para la gráfica con origen en (1, 0)
  vector_grafica <- c(1, 0, alineacion, distancia_normalizada)

  # Dibujar los vectores ideal y self
  plot(c(0, vector_ideal[1], vector_self[1]), c(0, vector_ideal[2], vector_self[2]), type = "n", xlab = "Coordenada X", ylab = "Coordenada Y", asp = 1, xlim = c(0, max(vector_ideal[1], vector_self[1]) + 1), ylim = c(0, max(vector_ideal[2], vector_self[2]) + 1))
  arrows(0, 0, vector_ideal[1], vector_ideal[2], col = "green", length = 0.1)
  arrows(0, 0, vector_self[1], vector_self[2], col = "orange", length = 0.1)
  text(vector_ideal[1] + 0.1, vector_ideal[2], "Vector Ideal", col = "green")
  text(vector_self[1] + 0.1, vector_self[2], "Vector Self", col = "orange")

  # Dibujar el único vector con origen en (1, 0) y destino en (alineacion, distancia_normalizada)
  plot(c(0, 1.5), c(0, 1.5), type = "n", xlab = "Alineación", ylab = "Distancia Normalizada", asp = 1, xlim = c(0, 1.5), ylim = c(0, 1.5))
  arrows(1, 0, alineacion, distancia_normalizada, col = "blue", length = 0.1)
  text(1.2, 0.1, "Vector (1,0) a (alineacion, distancia normalizada)", col = "blue")

  # Devolver los resultados
  return(list(distancia_normalizada = distancia_normalizada, alineacion = alineacion, grafica = vector_grafica))
}
