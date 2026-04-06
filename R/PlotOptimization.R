## PLOT OPTIMIZATION FUNCTIONS ##

# Enhanced label positioning algorithm inspired by ggrepel force-directed layout
.smart_label_positions <- function(x_coords, y_coords, labels = NULL,
                                   distance = 8, text_size = 11) {

  n <- length(x_coords)
  if (n == 0) return(list())

  # Remove any NA coordinates
  valid_idx <- !is.na(x_coords) & !is.na(y_coords)
  if (sum(valid_idx) == 0) return(list())

  x_coords <- x_coords[valid_idx]
  y_coords <- y_coords[valid_idx]
  if (!is.null(labels)) labels <- labels[valid_idx]
  n <- length(x_coords)

  # Initialize result data frame
  result <- data.frame(
    x = x_coords,
    y = y_coords,
    label = if(is.null(labels)) paste("Label", 1:n) else labels,
    xanchor = "center",
    yanchor = "middle",
    xshift = 0,
    yshift = 0,
    stringsAsFactors = FALSE
  )

  # Enhanced text dimension estimation with better precision
  char_width <- 0.007   # More conservative character width (was 0.006)
  char_height <- 0.013  # More conservative character height (was 0.012)

  # Calculate text box dimensions with character analysis
  text_widths <- sapply(result$label, function(txt) {
    char_count <- nchar(as.character(txt))
    # More sophisticated width calculation
    wide_chars <- length(gregexpr("[WMQwm@#%]", txt)[[1]])
    narrow_chars <- length(gregexpr("[ijlI.,;:]", txt)[[1]])
    base_width <- char_count * char_width * text_size / 11
    # Adjust for character width variations
    width_adjustment <- (wide_chars * 0.3) - (narrow_chars * 0.2)
    max(base_width + width_adjustment, base_width * 0.8)
  })
  text_heights <- rep(char_height * text_size / 11, n)

  # Adaptive parameters based on data density and plot constraints
  plot_margin <- 0.08  # 8% margin from edges (stricter boundary)

  # Dynamic distance adjustment based on density and overlaps -
  # more conservative
  point_density <- n / (diff(range(x_coords)) * diff(range(y_coords)))
  density_factor <- min(1.2, 0.8 / sqrt(point_density + 0.1))  # Less expansion

  base_distance <- distance * 1.5 * density_factor # Increased base distance multiplier
  min_distance <- max(2.5, base_distance * 0.1)  # Increased min distance
  max_distance <- base_distance * 0.6  # Increased max distance

  # Enhanced positioning with iterative force-based improvement (ggrepel-style)
  max_iterations <- min(50, n * 4)  # Increased iterations
  convergence_threshold <- 0.98  # Higher threshold (98% good positions)

  # System to lock excellent positions and prevent degradation
  locked_positions <- rep(FALSE, n)

  # Force-based parameters for spring system
  repulsion_force <- 0.03  # Increased repulsion (was 0.02)
  attraction_force <- 0.01
  damping <- 0.8  # Reduced damping for more movement (was 0.85)

  for (iteration in 1:max_iterations) {
    improved_count <- 0
    total_penalties <- 0

    # Detect dense clusters and apply zone-based optimization
    dense_zones <- .detect_label_clusters(
      result, text_widths, text_heights, x_coords, y_coords
    )

    # Smart ordering: process labels with fewer options first
    # (closer to edges/corners)
    label_order <- order(sapply(1:n, function(idx) {
      # Count how many directions are blocked by plot boundaries
      x <- x_coords[idx]
      y <- y_coords[idx]
      margin_penalty <- 0
      if (x < 0.3) margin_penalty <- margin_penalty + 1  # Left side
      if (x > 0.7) margin_penalty <- margin_penalty + 1  # Right side
      if (y < 0.3) margin_penalty <- margin_penalty + 1  # Bottom
      if (y > 0.7) margin_penalty <- margin_penalty + 1  # Top
      -margin_penalty
    }))

    for(i in label_order) {
      # Skip labels with locked excellent positions
      if(locked_positions[i]) next

      current_penalty <- .calculate_comprehensive_penalty(
        result, i, text_widths, text_heights, x_coords, y_coords,
        plot_margin
      )

      best_pos <- NULL
      min_penalty <- current_penalty
      found_excellent_orthogonal <- FALSE

      # Multi-strategy positioning with IMPROVED ORTHOGONAL PRIORITY
      # First, evaluate which orthogonal directions are available
      orthogonal_info <- .evaluate_orthogonal_spaces(
        x_coords[i], y_coords[i], result, i, text_widths, text_heights,
        x_coords, y_coords, plot_margin
      )

      # Apply cluster-aware positioning strategy
      in_dense_zone <- any(
        sapply(dense_zones, function(zone) i %in% zone)
      )

      if (in_dense_zone) {
        # For dense areas: prioritize the best available orthogonal positions
        strategies <- list(
          # ABSOLUTE priority for the specific best angles identified
          list(
            angles = orthogonal_info$best_angles[
              1:min(1, length(orthogonal_info$best_angles))
            ],
            dist_mult = 0.15, weight = 0.01
          ),
          # Secondary priority for all best angles
          list(
            angles = orthogonal_info$best_angles,
            dist_mult = 0.2, weight = 0.05
          ),
          # All orthogonal positions as backup
          list(
            angles = orthogonal_info$all_angles,
            dist_mult = 0.25, weight = 0.4
          ),
          # Mixed positions with orthogonal priority
          list(
            angles = c(
              orthogonal_info$all_angles, pi / 4, 3 * pi / 4,
              5 * pi / 4, 7 * pi / 4
            ),
            dist_mult = 0.3, weight = 0.8
          ),
          # Farther positions if needed
          list(
            angles = seq(0, 2 * pi, length.out = 16)[-16],
            dist_mult = 0.6, weight = 1.2
          )
        )
      } else {
        # For sparse areas: even stronger orthogonal preference
        strategies <- list(
          # ABSOLUTE priority for the specific best angle
          list(
            angles = orthogonal_info$best_angles[
              1:min(1, length(orthogonal_info$best_angles))
            ],
            dist_mult = 0.1, weight = 0.005
          ),
          # Secondary priority for all best angles
          list(
            angles = orthogonal_info$best_angles,
            dist_mult = 0.15, weight = 0.02
          ),
          # All orthogonal positions
          list(
            angles = orthogonal_info$all_angles,
            dist_mult = 0.2, weight = 0.3
          ),
          # Main diagonals as backup
          list(
            angles = c(pi / 4, 3 * pi / 4, 5 * pi / 4, 7 * pi / 4),
            dist_mult = 0.25, weight = 0.7
          ),
          # Last resort - any angle
          list(
            angles = seq(0, 2 * pi, length.out = 12)[-12],
            dist_mult = 0.4, weight = 1.1
          )
        )
      }

      for (strategy in strategies) {
        strategy_found_good <- FALSE

        for (angle in strategy$angles) {
          # Calculate position with current strategy
          test_distance <- min_distance +
            (max_distance - min_distance) * strategy$dist_mult

          # Check if this is an orthogonal angle for reduced jitter
          orthogonal_angles <- c(0, pi / 2, pi, 3 * pi / 2)
          is_orthogonal <- any(abs(angle - orthogonal_angles) < 0.1)

          # Apply minimal jitter for orthogonal angles, more for diagonals
          jitter_amount <- ifelse(is_orthogonal, 0.0005, 0.002)
          jitter_x <- rnorm(1, 0, jitter_amount)
          jitter_y <- rnorm(1, 0, jitter_amount)

          # Improved distance calculation - adaptive for edge positions
          base_dist <- ifelse(is_orthogonal, 0.03, 0.025) # Increased from 0.02, 0.015

          # Reduce distance for points near boundaries to avoid violations
          distance_from_edges <- min(
            x_coords[i] - plot_margin,
            (1 - plot_margin) - x_coords[i],
            y_coords[i] - plot_margin,
            (1 - plot_margin) - y_coords[i]
          )

          # Scale down distance if too close to edges
          if (distance_from_edges < 0.08) {
            edge_factor <- max(0.3, distance_from_edges / 0.08)
            base_dist <- base_dist * edge_factor
          }

          actual_distance <- base_dist + test_distance * 0.01

          label_x <- x_coords[i] + cos(angle) * actual_distance + jitter_x
          label_y <- y_coords[i] + sin(angle) * actual_distance + jitter_y

          # Estimate actual label boundaries based on text size and anchor
          estimated_label_width <- text_widths[i]
          estimated_label_height <- text_heights[i]

          # Calculate label bounds for boundary checking
          label_left <- label_x - estimated_label_width / 2
          label_right <- label_x + estimated_label_width / 2
          label_bottom <- label_y - estimated_label_height / 2
          label_top <- label_y + estimated_label_height / 2

          # Strict boundary enforcement - reject if ANY part goes outside
          if (label_left < plot_margin || label_right > (1 - plot_margin) ||
              label_bottom < plot_margin || label_top > (1 - plot_margin)) {
            next  # Skip this position completely
          }

          # Additional safety check - ensure center is well within bounds
          if (label_x < plot_margin + estimated_label_width / 2 ||
              label_x > (1 - plot_margin - estimated_label_width / 2) ||
              label_y < plot_margin + estimated_label_height / 2 ||
              label_y > (1 - plot_margin - estimated_label_height / 2)) {
            next  # Skip this position too
          }

          # Create temporary result for testing
          test_result <- result
          test_result <- .update_label_position_precise(
            test_result, i, label_x, label_y, x_coords[i], y_coords[i],
            test_distance, text_widths[i], text_heights[i]
          )

          # Calculate comprehensive penalty
          penalty <- .calculate_comprehensive_penalty(
            test_result, i, text_widths, text_heights, x_coords, y_coords,
            plot_margin
          ) * strategy$weight  # Weight by strategy preference

          # Enhanced bonus system for orthogonal positions
          if (is_orthogonal) {
            # Check if this angle is in the "best" available angles
            is_best_angle <- angle %in% orthogonal_info$best_angles
            # Check if this is the SINGLE best angle (first in the list)
            is_top_choice <- length(orthogonal_info$best_angles) > 0 &&
              abs(angle - orthogonal_info$best_angles[1]) < 0.1

            if (is_top_choice && penalty < 5) {
              penalty <- penalty * 0.02  # 98% bonus for the absolute best choice
            } else if (is_best_angle && penalty < 5) {
              penalty <- penalty * 0.05  # 95% bonus for excellent orthogonal
            } else if (penalty < 8) {
              penalty <- penalty * 0.2  # 80% bonus for good orthogonal
            } else if (penalty < 15) {
              penalty <- penalty * 0.5  # 50% bonus for acceptable orthogonal
            }
          }

          # Penalty for distance from point (prefer closer labels)
          distance_penalty <- actual_distance * 8  # Reduced distance penalty
          penalty <- penalty + distance_penalty

          if (penalty < min_penalty) {
            min_penalty <- penalty
            best_pos <- list(
              x = label_x, y = label_y, angle = angle,
              strategy = strategy, distance = test_distance,
              is_orthogonal = is_orthogonal
            )
            strategy_found_good <- TRUE
            improved_count <- improved_count + 1

            # If found a good orthogonal position, use it immediately
            if (is_orthogonal && penalty < 1.0) {
              found_excellent_orthogonal <- TRUE
              break  # Exit angle loop immediately
            }
          }
        }

        # Early exit if found excellent orthogonal position
        if (found_excellent_orthogonal) {
          break  # Exit strategy loop immediately
        }

        # Early exit if found good orthogonal position
        if (strategy_found_good && !is.null(best_pos) &&
            best_pos$is_orthogonal && min_penalty < 3.0) {
          break  # Exit strategy loop immediately
        }

        # Early exit if close strategy found good position
        if (strategy_found_good && strategy$dist_mult <= 0.2) break
      }

      # Apply best position if found
      if (!is.null(best_pos)) {
        result <- .update_label_position_precise(
          result, i, best_pos$x, best_pos$y,
          x_coords[i], y_coords[i], best_pos$distance,
          text_widths[i], text_heights[i]
        )

        # Lock positions that are excellent and orthogonal
        if (best_pos$is_orthogonal && min_penalty < 1.0) {
          locked_positions[i] <- TRUE
          improved_count <- improved_count + 1
        }
      }

      total_penalties <- total_penalties + min_penalty
    }

    # Convergence check - stop if most positions are good
    improvement_ratio <- improved_count / n
    if (improvement_ratio >= convergence_threshold) break

    # Adaptive parameter adjustment for next iteration
    if (iteration > 5 && improvement_ratio < 0.3) {
      # If not improving much, slightly relax distance constraints
      max_distance <- min(max_distance * 1.1, base_distance)
    }
  }

  # Final boundary enforcement and simple overlap resolution
  result <- .force_labels_within_bounds(
    result, text_widths, text_heights, x_coords, y_coords, plot_margin
  )
  result <- .final_overlap_cleanup(
    result, text_widths, text_heights, x_coords, y_coords, plot_margin
  )

  result
}

# Final cleanup to ensure no overlaps remain
.final_overlap_cleanup <- function(
    result, text_widths, text_heights, x_coords, y_coords, plot_margin) {
  n <- nrow(result)

  for (attempt in 1:3) {
    any_overlap <- FALSE

    for (i in 1:n) {
      current_x <- result$x[i] + result$xshift[i] * 0.0008
      current_y <- result$y[i] + result$yshift[i] * 0.0008

      # Check for overlaps with other labels
      for (j in (i + 1):n) {
        if (j > n) break

        other_x <- result$x[j] + result$xshift[j] * 0.0008
        other_y <- result$y[j] + result$yshift[j] * 0.0008

        # Check overlap
        x_overlap <- min(
          current_x + text_widths[i] / 2, other_x + text_widths[j] / 2
        ) -
          max(current_x - text_widths[i] / 2, other_x - text_widths[j] / 2)
        y_overlap <- min(
          current_y + text_heights[i] / 2, other_y + text_heights[j] / 2
        ) -
          max(current_y - text_heights[i] / 2, other_y - text_heights[j] / 2)

        if (x_overlap > 0.001 && y_overlap > 0.001) {
          any_overlap <- TRUE

          # Move label j to a nearby non-overlapping position
          for (angle in seq(0, 2 * pi, length.out = 8)) {
            for (dist in c(0.02, 0.03)) {
              new_x <- x_coords[j] + cos(angle) * dist
              new_y <- y_coords[j] + sin(angle) * dist

              # Check bounds
              if (new_x - text_widths[j] / 2 >= plot_margin &&
                  new_x + text_widths[j] / 2 <= (1 - plot_margin) &&
                  new_y - text_heights[j] / 2 >= plot_margin &&
                  new_y + text_heights[j] / 2 <= (1 - plot_margin)) {

                # Check no overlap with current label
                test_x_overlap <- min(
                  current_x + text_widths[i] / 2,
                  new_x + text_widths[j] / 2
                ) -
                  max(
                    current_x - text_widths[i] / 2,
                    new_x - text_widths[j] / 2
                  )
                test_y_overlap <- min(
                  current_y + text_heights[i] / 2,
                  new_y + text_heights[j] / 2
                ) -
                  max(
                    current_y - text_heights[i] / 2,
                    new_y - text_heights[j] / 2
                  )

                if (test_x_overlap <= 0.001 || test_y_overlap <= 0.001) {
                  # Good position - update
                  result <- .update_label_position_precise(
                    result, j, new_x, new_y, x_coords[j], y_coords[j],
                    dist * 100, text_widths[j], text_heights[j]
                  )
                  break
                }
              }
            }
          }
        }
      }
    }

    if (!any_overlap) break
  }

  result
}

# Detect clusters of labels in dense areas for targeted optimization
.detect_label_clusters <- function(
    result, text_widths, text_heights, x_coords, y_coords) {
  n <- nrow(result)
  clusters <- list()
  visited <- rep(FALSE, n)

  for (i in 1:n) {
    if (visited[i]) next

    # Start a new cluster
    cluster <- c(i)
    visited[i] <- TRUE

    current_x <- result$x[i] + result$xshift[i] * 0.0008
    current_y <- result$y[i] + result$yshift[i] * 0.0008

    # Find nearby labels to add to cluster
    for (j in 1:n) {
      if (visited[j] || i == j) next

      other_x <- result$x[j] + result$xshift[j] * 0.0008
      other_y <- result$y[j] + result$yshift[j] * 0.0008

      # Distance threshold for clustering
      dist_threshold <- max(
        text_widths[i], text_heights[i],
        text_widths[j], text_heights[j]
      ) * 2.5
      distance <- sqrt((current_x - other_x)^2 + (current_y - other_y)^2)

      if (distance < dist_threshold) {
        cluster <- c(cluster, j)
        visited[j] <- TRUE
      }
    }

    # Only consider it a cluster if it has multiple labels
    if (length(cluster) > 1) {
      clusters[[length(clusters) + 1]] <- cluster
    }
  }

  clusters
}

# Apply spring forces for final optimization (ggrepel-inspired)
.apply_spring_forces <- function(
    result, text_widths, text_heights, x_coords, y_coords, plot_margin,
    attraction_force, repulsion_force, damping) {
  n <- nrow(result)
  max_force_iterations <- 15

  # Initialize velocities
  velocities_x <- rep(0, n)
  velocities_y <- rep(0, n)

  for (iteration in 1:max_force_iterations) {
    forces_x <- rep(0, n)
    forces_y <- rep(0, n)

    # Calculate forces for each label
    for (i in 1:n) {
      current_x <- result$x[i] + result$xshift[i] * 0.0008
      current_y <- result$y[i] + result$yshift[i] * 0.0008

      # Attraction force to original marker position (spring force)
      dx_to_marker <- x_coords[i] - current_x
      dy_to_marker <- y_coords[i] - current_y
      marker_distance <- sqrt(dx_to_marker^2 + dy_to_marker^2)

      if (marker_distance > 0.02) {  # Only if not too close
        forces_x[i] <- forces_x[i] + dx_to_marker * attraction_force
        forces_y[i] <- forces_y[i] + dy_to_marker * attraction_force
      }

      # Repulsion forces from other labels
      for (j in 1:n) {
        if (i == j) next

        other_x <- result$x[j] + result$xshift[j] * 0.0008
        other_y <- result$y[j] + result$yshift[j] * 0.0008

        dx <- current_x - other_x
        dy <- current_y - other_y
        distance <- sqrt(dx^2 + dy^2)

        # Minimum safe distance based on label sizes
        min_distance <- (
          max(text_widths[i], text_heights[i]) +
            max(text_widths[j], text_heights[j])
        ) * 0.8

        if (distance < min_distance && distance > 0.001) {
          # Repulsion force (inverse square law)
          force_magnitude <- repulsion_force / (distance^2 + 0.001)
          force_magnitude <- min(force_magnitude, 0.01)  # Cap maximum force

          forces_x[i] <- forces_x[i] + (dx / distance) * force_magnitude
          forces_y[i] <- forces_y[i] + (dy / distance) * force_magnitude
        }
      }

      # Boundary repulsion forces
      boundary_force_strength <- 0.005
      margin_distance <- 0.05

      if(current_x < plot_margin + margin_distance) {
        forces_x[i] <- forces_x[i] + boundary_force_strength / (current_x - plot_margin + 0.001)
      }
      if(current_x > (1 - plot_margin - margin_distance)) {
        forces_x[i] <- forces_x[i] - boundary_force_strength / ((1 - plot_margin) - current_x + 0.001)
      }
      if(current_y < plot_margin + margin_distance) {
        forces_y[i] <- forces_y[i] + boundary_force_strength / (current_y - plot_margin + 0.001)
      }
      if(current_y > (1 - plot_margin - margin_distance)) {
        forces_y[i] <- forces_y[i] - boundary_force_strength / ((1 - plot_margin) - current_y + 0.001)
      }
    }

    # Update velocities and positions
    for(i in 1:n) {
      # Update velocity with damping
      velocities_x[i] <- (velocities_x[i] + forces_x[i]) * damping
      velocities_y[i] <- (velocities_y[i] + forces_y[i]) * damping

      # Update position
      new_x <- result$x[i] + result$xshift[i] * 0.0008 + velocities_x[i]
      new_y <- result$y[i] + result$yshift[i] * 0.0008 + velocities_y[i]

      # Ensure within boundaries
      new_x <- max(plot_margin + text_widths[i]/2, min((1 - plot_margin) - text_widths[i]/2, new_x))
      new_y <- max(plot_margin + text_heights[i]/2, min((1 - plot_margin) - text_heights[i]/2, new_y))

      # Update result
      result <- .update_label_position_precise(
        result, i, new_x, new_y, x_coords[i], y_coords[i], 8, text_widths[i], text_heights[i]
      )
    }
  }

  return(result)
}

# Resolve any remaining collisions after main algorithm
.resolve_remaining_collisions <- function(result, text_widths, text_heights, x_coords, y_coords, plot_margin) {
  n <- nrow(result)
  max_collision_resolution_attempts <- 10

  for(attempt in 1:max_collision_resolution_attempts) {
    collisions_found <- FALSE

    for(i in 1:n) {
      # Check if this label collides with any other
      current_x <- result$x[i] + result$xshift[i] * 0.0008
      current_y <- result$y[i] + result$yshift[i] * 0.0008
      current_w <- text_widths[i]
      current_h <- text_heights[i]

      for(j in (i+1):n) {
        if(j > n) break

        other_x <- result$x[j] + result$xshift[j] * 0.0008
        other_y <- result$y[j] + result$yshift[j] * 0.0008
        other_w <- text_widths[j]
        other_h <- text_heights[j]

        # Check for collision with expanded safety margin
        margin <- 0.008
        x_overlap <- max(0, min(current_x + current_w/2 + margin, other_x + other_w/2 + margin) -
                              max(current_x - current_w/2 - margin, other_x - other_w/2 - margin))
        y_overlap <- max(0, min(current_y + current_h/2 + margin, other_y + other_h/2 + margin) -
                              max(current_y - current_h/2 - margin, other_y - other_h/2 - margin))

        if(x_overlap > 0 && y_overlap > 0) {
          collisions_found <- TRUE

          # Move the label with higher index (j) away from the collision
          # Calculate direction to move away
          dx <- other_x - current_x
          dy <- other_y - current_y

          if(abs(dx) > abs(dy)) {
            # Move horizontally
            move_distance <- (current_w + other_w)/2 + margin + 0.005
            new_x <- current_x + sign(dx) * move_distance
          } else {
            # Move vertically
            move_distance <- (current_h + other_h)/2 + margin + 0.005
            new_y <- current_y + sign(dy) * move_distance
            new_x <- other_x
          }

          # Ensure the new position is within bounds
          if(abs(dx) > abs(dy)) {
            if(new_x - other_w/2 >= plot_margin && new_x + other_w/2 <= (1-plot_margin)) {
              result <- .update_label_position_precise(
                result, j, new_x, other_y, x_coords[j], y_coords[j], 8, other_w, other_h
              )
            }
          } else {
            if(new_y - other_h/2 >= plot_margin && new_y + other_h/2 <= (1-plot_margin)) {
              result <- .update_label_position_precise(
                result, j, new_x, new_y, x_coords[j], y_coords[j], 8, other_w, other_h
              )
            }
          }
        }
      }
    }

    # If no more collisions found, we're done
    if(!collisions_found) break
  }

  return(result)
}

# Force labels within plot boundaries as final safety measure
.force_labels_within_bounds <- function(result, text_widths, text_heights, x_coords, y_coords, plot_margin) {
  n <- nrow(result)

  for(i in 1:n) {
    # Calculate current estimated position
    current_x <- result$x[i] + result$xshift[i] * 0.0008
    current_y <- result$y[i] + result$yshift[i] * 0.0008
    label_w <- text_widths[i]
    label_h <- text_heights[i]

    # Check if label extends outside boundaries
    left_bound <- current_x - label_w/2
    right_bound <- current_x + label_w/2
    bottom_bound <- current_y - label_h/2
    top_bound <- current_y + label_h/2

    # Flag if any part is outside
    outside <- FALSE
    new_x <- current_x
    new_y <- current_y

    # Clamp to boundaries with margin
    if(left_bound < plot_margin) {
      new_x <- plot_margin + label_w/2
      outside <- TRUE
    } else if(right_bound > (1 - plot_margin)) {
      new_x <- (1 - plot_margin) - label_w/2
      outside <- TRUE
    }

    if(bottom_bound < plot_margin) {
      new_y <- plot_margin + label_h/2
      outside <- TRUE
    } else if(top_bound > (1 - plot_margin)) {
      new_y <- (1 - plot_margin) - label_h/2
      outside <- TRUE
    }

    # If we had to move the label, update with corrected position
    if(outside) {
      # Find the best anchor/shift combination for the corrected position
      result <- .update_label_position_precise(
        result, i, new_x, new_y, x_coords[i], y_coords[i],
        8, label_w, label_h  # Use default distance
      )
    }
  }

  return(result)
}

# Simplified penalty calculation focused on eliminating overlaps while staying close
.calculate_overlap_focused_penalty <- function(result, index, text_widths, text_heights, x_coords, y_coords, plot_margin = 0.08) {
  penalty <- 0

  # Current label position estimation
  current_x <- result$x[index] + result$xshift[index] * 0.0008
  current_y <- result$y[index] + result$yshift[index] * 0.0008
  current_w <- text_widths[index]
  current_h <- text_heights[index]

  n <- nrow(result)

  # Focus on label-label overlaps with zero tolerance
  for(j in 1:n) {
    if(j == index) next

    other_x <- result$x[j] + result$xshift[j] * 0.0008
    other_y <- result$y[j] + result$yshift[j] * 0.0008
    other_w <- text_widths[j]
    other_h <- text_heights[j]

    # Zero tolerance overlap detection with small buffer
    buffer <- 0.003  # Very small buffer for visual separation
    x_overlap <- max(0, min(current_x + current_w/2, other_x + other_w/2) -
                          max(current_x - current_w/2, other_x - other_w/2) + buffer)
    y_overlap <- max(0, min(current_y + current_h/2, other_y + other_h/2) -
                          max(current_y - current_h/2, other_y - other_h/2) + buffer)

    # Massive penalty for ANY overlap
    if(x_overlap > 0 && y_overlap > 0) {
      penalty <- penalty + 50000  # Fixed massive penalty
    }

    # Proximity penalty to maintain spacing
    label_dist <- sqrt((current_x - other_x)^2 + (current_y - other_y)^2)
    min_spacing <- (max(current_w, current_h) + max(other_w, other_h)) * 0.6
    if(label_dist < min_spacing) {
      penalty <- penalty + (min_spacing - label_dist) * 5000
    }
  }

  # Marker collision check - moderate penalty for other markers
  for(j in seq_along(x_coords)) {
    marker_dist <- sqrt((current_x - x_coords[j])^2 + (current_y - y_coords[j])^2)

    if(j != index) {
      # Other markers - avoid collision
      min_marker_distance <- max(current_w, current_h) * 0.5 + 0.012
      if(marker_dist < min_marker_distance) {
        penalty <- penalty + (min_marker_distance - marker_dist) * 3000
      }
    } else {
      # Own marker - gentle preference for proximity
      if(marker_dist > 0.05) {
        penalty <- penalty + (marker_dist - 0.05) * 100
      }
    }
  }

  return(penalty)
}

# Simplified penalty calculation focused on eliminating overlaps while staying close
.calculate_overlap_focused_penalty <- function(result, index, text_widths, text_heights, x_coords, y_coords, plot_margin = 0.08) {
  penalty <- 0

  # Current label position estimation
  current_x <- result$x[index] + result$xshift[index] * 0.0008
  current_y <- result$y[index] + result$yshift[index] * 0.0008
  current_w <- text_widths[index]
  current_h <- text_heights[index]

  n <- nrow(result)

  # Focus on label-label overlaps with zero tolerance
  for(j in 1:n) {
    if(j == index) next

    other_x <- result$x[j] + result$xshift[j] * 0.0008
    other_y <- result$y[j] + result$yshift[j] * 0.0008
    other_w <- text_widths[j]
    other_h <- text_heights[j]

    # Zero tolerance overlap detection with small buffer
    buffer <- 0.003  # Very small buffer for visual separation
    x_overlap <- max(0, min(current_x + current_w/2, other_x + other_w/2) -
                          max(current_x - current_w/2, other_x - other_w/2) + buffer)
    y_overlap <- max(0, min(current_y + current_h/2, other_y + other_h/2) -
                          max(current_y - current_h/2, other_y - other_h/2) + buffer)

    # Massive penalty for ANY overlap
    if(x_overlap > 0 && y_overlap > 0) {
      penalty <- penalty + 50000  # Fixed massive penalty
    }

    # Proximity penalty to maintain spacing
    label_dist <- sqrt((current_x - other_x)^2 + (current_y - other_y)^2)
    min_spacing <- (max(current_w, current_h) + max(other_w, other_h)) * 0.6
    if(label_dist < min_spacing) {
      penalty <- penalty + (min_spacing - label_dist) * 5000
    }
  }

  # Marker collision check - moderate penalty for other markers
  for(j in seq_along(x_coords)) {
    marker_dist <- sqrt((current_x - x_coords[j])^2 + (current_y - y_coords[j])^2)

    if(j != index) {
      # Other markers - avoid collision
      min_marker_distance <- max(current_w, current_h) * 0.5 + 0.012
      if(marker_dist < min_marker_distance) {
        penalty <- penalty + (min_marker_distance - marker_dist) * 3000
      }
    } else {
      # Own marker - gentle preference for proximity
      if(marker_dist > 0.05) {
        penalty <- penalty + (marker_dist - 0.05) * 100
      }
    }
  }

  return(penalty)
}

# Simple but effective overlap resolution
.resolve_overlaps_simple <- function(result, text_widths, text_heights, x_coords, y_coords, plot_margin) {
  n <- nrow(result)
  max_attempts <- 5

  for(attempt in 1:max_attempts) {
    overlaps_found <- FALSE

    for(i in 1:n) {
      current_x <- result$x[i] + result$xshift[i] * 0.0008
      current_y <- result$y[i] + result$yshift[i] * 0.0008
      current_w <- text_widths[i]
      current_h <- text_heights[i]

      for(j in (i+1):n) {
        if(j > n) break

        other_x <- result$x[j] + result$xshift[j] * 0.0008
        other_y <- result$y[j] + result$yshift[j] * 0.0008
        other_w <- text_widths[j]
        other_h <- text_heights[j]

        # Check for overlap
        x_overlap <- max(0, min(current_x + current_w/2, other_x + other_w/2) -
                              max(current_x - current_w/2, other_x - other_w/2))
        y_overlap <- max(0, min(current_y + current_h/2, other_y + other_h/2) -
                              max(current_y - current_h/2, other_y - other_h/2))

        if(x_overlap > 0 && y_overlap > 0) {
          overlaps_found <- TRUE

          # Move the second label (j) away - find closest valid position
          best_new_pos <- NULL
          min_distance_to_marker <- Inf

          # Try positions around its marker at small distance
          for(angle in seq(0, 2*pi, length.out = 12)[-12]) {
            for(dist in c(0.02, 0.03, 0.04)) {
              new_x <- x_coords[j] + cos(angle) * dist
              new_y <- y_coords[j] + sin(angle) * dist

              # Check if within bounds
              if(new_x - other_w/2 < plot_margin || new_x + other_w/2 > (1-plot_margin) ||
                 new_y - other_h/2 < plot_margin || new_y + other_h/2 > (1-plot_margin)) {
                next
              }

              # Check if this position avoids overlap with current label
              test_x_overlap <- max(0, min(current_x + current_w/2, new_x + other_w/2) -
                                         max(current_x - current_w/2, new_x - other_w/2))
              test_y_overlap <- max(0, min(current_y + current_h/2, new_y + other_h/2) -
                                         max(current_y - current_h/2, new_y - other_h/2))

              if(test_x_overlap <= 0.001 || test_y_overlap <= 0.001) {
                # No overlap - check distance to marker
                dist_to_marker <- sqrt((new_x - x_coords[j])^2 + (new_y - y_coords[j])^2)
                if(dist_to_marker < min_distance_to_marker) {
                  min_distance_to_marker <- dist_to_marker
                  best_new_pos <- list(x = new_x, y = new_y)
                }
              }
            }
          }

          # Apply the best position if found
          if(!is.null(best_new_pos)) {
            result <- .update_label_position_precise(
              result, j, best_new_pos$x, best_new_pos$y, x_coords[j], y_coords[j],
              min_distance_to_marker * 100, other_w, other_h
            )
          }
        }
      }
    }

    if(!overlaps_found) break
  }

  return(result)
}

# Enhanced helper functions inspired by ggrepel algorithms ---------------------

# Precise label position update with boundary awareness
.update_label_position_precise <- function(result, index, label_x, label_y, point_x, point_y, distance, label_width, label_height) {
  # Calculate offset from marker center to label position
  dx <- label_x - point_x
  dy <- label_y - point_y

  # Enhanced anchor system with tighter control
  abs_dx <- abs(dx)
  abs_dy <- abs(dy)

  # Determine primary direction with higher threshold for precision
  if(abs_dx > abs_dy * 1.5) {
    # Horizontal positioning - label to left/right of point
    result$xanchor[index] <- if(dx > 0) "left" else "right"
    result$yanchor[index] <- "middle"
    # Tighter spacing - stay close to marker
    result$xshift[index] <- if(dx > 0) max(4, distance * 0.4) else -max(4, distance * 0.4)
    result$yshift[index] <- dy * 200  # Minimal vertical adjustment
  } else if(abs_dy > abs_dx * 1.5) {
    # Vertical positioning - label above/below point
    result$xanchor[index] <- "center"
    result$yanchor[index] <- if(dy > 0) "bottom" else "top"
    result$xshift[index] <- dx * 200  # Minimal horizontal adjustment
    result$yshift[index] <- if(dy > 0) max(4, distance * 0.4) else -max(4, distance * 0.4)
  } else {
    # Diagonal positioning - keep it tight
    result$xanchor[index] <- if(dx > 0) "left" else "right"
    result$yanchor[index] <- if(dy > 0) "bottom" else "top"
    shift_amount <- max(3, distance * 0.3)
    result$xshift[index] <- if(dx > 0) shift_amount else -shift_amount
    result$yshift[index] <- if(dy > 0) shift_amount else -shift_amount
  }

  return(result)
}

# Backwards compatibility function
.update_label_position <- function(result, index, label_x, label_y, point_x, point_y, distance) {
  # Use default label dimensions for backwards compatibility
  default_width <- 0.05
  default_height <- 0.02
  return(.update_label_position_precise(result, index, label_x, label_y, point_x, point_y, distance, default_width, default_height))
}

# Comprehensive penalty calculation with enhanced collision detection
.calculate_comprehensive_penalty <- function(result, index, text_widths, text_heights, x_coords, y_coords, plot_margin = 0.08) {
  penalty <- 0

  # Current label position estimation with improved accuracy
  current_x <- result$x[index] + result$xshift[index] * 0.0008
  current_y <- result$y[index] + result$yshift[index] * 0.0008
  current_w <- text_widths[index]
  current_h <- text_heights[index]

  # Strict boundary checking with margin
  if(current_x - current_w/2 < plot_margin) penalty <- penalty + 8000
  if(current_x + current_w/2 > (1 - plot_margin)) penalty <- penalty + 8000
  if(current_y - current_h/2 < plot_margin) penalty <- penalty + 8000
  if(current_y + current_h/2 > (1 - plot_margin)) penalty <- penalty + 8000

  n <- nrow(result)

  # Enhanced label-label collision detection with larger safety margins
  for(j in 1:n) {
    if(j == index) next

    # Other label position with improved estimation
    other_x <- result$x[j] + result$xshift[j] * 0.0008
    other_y <- result$y[j] + result$yshift[j] * 0.0008
    other_w <- text_widths[j]
    other_h <- text_heights[j]

    # Increased safety margins to prevent any visual overlap
    safety_margin_x <- max(current_w, other_w) * 0.15  # 15% of larger width
    safety_margin_y <- max(current_h, other_h) * 0.15  # 15% of larger height

    # Enhanced box overlap calculation with safety margins
    x_overlap <- max(0, min(current_x + current_w/2 + safety_margin_x, other_x + other_w/2 + safety_margin_x) -
                          max(current_x - current_w/2 - safety_margin_x, other_x - other_w/2 - safety_margin_x))
    y_overlap <- max(0, min(current_y + current_h/2 + safety_margin_y, other_y + other_h/2 + safety_margin_y) -
                          max(current_y - current_h/2 - safety_margin_y, other_y - other_h/2 - safety_margin_y))

    if(x_overlap > 0 && y_overlap > 0) {
      # Massive penalty for any overlap - make this prohibitively expensive
      overlap_area <- x_overlap * y_overlap
      penalty <- penalty + overlap_area * 50000 * (4.0 ^ (overlap_area * 600))
    }

    # Progressive proximity penalty with larger minimum distances
    label_dist <- sqrt((current_x - other_x)^2 + (current_y - other_y)^2)
    min_safe_distance <- (max(current_w, current_h) + max(other_w, other_h)) * 0.8  # Increased from 0.6
    if(label_dist < min_safe_distance) {
      proximity_penalty <- (min_safe_distance - label_dist)^3 * 2000  # Increased penalty
      penalty <- penalty + proximity_penalty
    }
  }

  # Enhanced marker collision check for ALL markers
  for(j in seq_along(x_coords)) {
    marker_dist <- sqrt((current_x - x_coords[j])^2 + (current_y - y_coords[j])^2)

    # Dynamic minimum distance based on label size
    min_marker_distance <- max(current_w, current_h) * 0.8 + 0.025  # Increased buffer

    if(marker_dist < min_marker_distance) {
      if(j != index) {  # Other markers
        collision_penalty <- (min_marker_distance - marker_dist)^2 * 2500
      } else {  # Own marker - allow closer proximity but not too close
        own_marker_min <- max(current_w, current_h) * 0.3 + 0.010
        if(marker_dist < own_marker_min) {
          collision_penalty <- (own_marker_min - marker_dist)^2 * 800
        } else {
          collision_penalty <- 0
        }
      }
      penalty <- penalty + collision_penalty
    }

    # Distance preference - prefer staying reasonably close to own marker
    if(j == index && marker_dist > 0.08) {  # Slightly increased max distance
      distance_penalty <- (marker_dist - 0.08)^2 * 60
      penalty <- penalty + distance_penalty
    }
  }

  return(penalty)
}

# Backwards compatibility function
.calculate_enhanced_penalty <- function(result, index, text_widths, text_heights, x_coords, y_coords) {
  return(.calculate_comprehensive_penalty(result, index, text_widths, text_heights, x_coords, y_coords, 0.05))
}

# Boundary force calculation (inspired by ggrepel's put_within_bounds)
.calculate_boundary_penalty <- function(label_x, label_y, force_strength) {
  penalty <- 0
  margin <- 0.05  # 5% margin from edges

  # Boundary repulsion forces (exponential near edges)
  if(label_x < margin) {
    penalty <- penalty + (margin - label_x)^2 * 1000 / force_strength
  }
  if(label_x > 1 - margin) {
    penalty <- penalty + (label_x - (1 - margin))^2 * 1000 / force_strength
  }
  if(label_y < margin) {
    penalty <- penalty + (margin - label_y)^2 * 1000 / force_strength
  }
  if(label_y > 1 - margin) {
    penalty <- penalty + (label_y - (1 - margin))^2 * 1000 / force_strength
  }

  return(penalty)
}

# Diversity penalty to avoid position clustering
.calculate_diversity_penalty <- function(result, index, angle) {
  penalty <- 0
  n <- nrow(result)

  if(n <= 1 || index <= 1) return(0)

  # Calculate position type from angle
  current_type <- .enhanced_angle_to_position_type(angle)

  # Count usage of similar positions
  type_count <- 0
  for(j in 1:(index-1)) {
    if(j <= n) {  # Safety check
      other_type <- paste(result$xanchor[j], result$yanchor[j])
      if(length(other_type) > 0 && other_type == current_type) {
        type_count <- type_count + 1
      }
    }
  }

  # Penalty increases exponentially with overuse (like ggrepel's position balancing)
  if(type_count > 0) {
    penalty <- penalty + type_count^1.5 * 3
  }

  return(penalty)
}

# Enhanced position type classification
.enhanced_angle_to_position_type <- function(angle) {
  angle_deg <- (angle * 180 / pi) %% 360

  # 8-direction classification for better precision
  if(angle_deg >= 337.5 || angle_deg < 22.5) return("left middle")
  if(angle_deg >= 22.5 && angle_deg < 67.5) return("left bottom")
  if(angle_deg >= 67.5 && angle_deg < 112.5) return("center bottom")
  if(angle_deg >= 112.5 && angle_deg < 157.5) return("right bottom")
  if(angle_deg >= 157.5 && angle_deg < 202.5) return("right middle")
  if(angle_deg >= 202.5 && angle_deg < 247.5) return("right top")
  if(angle_deg >= 247.5 && angle_deg < 292.5) return("center top")
  if(angle_deg >= 292.5 && angle_deg < 337.5) return("left top")

  return("center middle")
}

# Evaluate available orthogonal spaces for intelligent positioning
# Enhanced evaluation of available orthogonal spaces for intelligent positioning
.evaluate_orthogonal_spaces <- function(point_x, point_y, result, current_idx,
                                        text_widths, text_heights, x_coords, y_coords,
                                        plot_margin) {

  # Define orthogonal angles and their space requirements
  orthogonal_angles <- c(0, pi/2, pi, 3*pi/2)  # Right, Top, Left, Bottom
  angle_names <- c("right", "top", "left", "bottom")

  # Evaluate each orthogonal direction
  space_scores <- numeric(4)
  names(space_scores) <- angle_names

  for(j in 1:length(orthogonal_angles)) {
    angle <- orthogonal_angles[j]
    angle_name <- angle_names[j]

    # Use multiple test distances to find optimal positioning
    min_penalty <- Inf
    best_distance <- 0.025

    for(test_distance in c(0.02, 0.025, 0.03, 0.035, 0.04)) {
      test_x <- point_x + cos(angle) * test_distance
      test_y <- point_y + sin(angle) * test_distance

      # Check boundary constraints with actual label dimensions
      label_width <- text_widths[current_idx]
      label_height <- text_heights[current_idx]

      # Calculate precise label bounds
      label_left <- test_x - label_width/2
      label_right <- test_x + label_width/2
      label_bottom <- test_y - label_height/2
      label_top <- test_y + label_height/2

          # Strict boundary checking with tolerance for extreme edge cases
          boundary_penalty <- 0
          margin_tolerance <- 0.01  # Allow small violations for extreme cases

          if(label_left < plot_margin - margin_tolerance) boundary_penalty <- boundary_penalty + 100
          if(label_right > (1 - plot_margin + margin_tolerance)) boundary_penalty <- boundary_penalty + 100
          if(label_bottom < plot_margin - margin_tolerance) boundary_penalty <- boundary_penalty + 100
          if(label_top > (1 - plot_margin + margin_tolerance)) boundary_penalty <- boundary_penalty + 100

          # For very minor violations, add smaller penalty instead of complete rejection
          if(boundary_penalty == 0) {
            minor_violation <- 0
            if(label_left < plot_margin) minor_violation <- minor_violation + (plot_margin - label_left) * 1000
            if(label_right > (1 - plot_margin)) minor_violation <- minor_violation + (label_right - (1 - plot_margin)) * 1000
            if(label_bottom < plot_margin) minor_violation <- minor_violation + (plot_margin - label_bottom) * 1000
            if(label_top > (1 - plot_margin)) minor_violation <- minor_violation + (label_top - (1 - plot_margin)) * 1000
            boundary_penalty <- minor_violation
          }

          # Skip only for major boundary violations
          if(boundary_penalty >= 100) next

          # Enhanced collision detection with actual geometric overlap
          collision_penalty <- 0

      for(k in 1:length(x_coords)) {
        if(k == current_idx) next

        # Distance to other points (markers)
        point_dist <- sqrt((test_x - x_coords[k])^2 + (test_y - y_coords[k])^2)

        # Strong penalty for being too close to other markers
        if(point_dist < 0.03) collision_penalty <- collision_penalty + 50
        else if(point_dist < 0.05) collision_penalty <- collision_penalty + 20

        # Distance and overlap to other labels (if they exist)
        if(k <= nrow(result)) {
          other_label_x <- x_coords[k] + result$xshift[k] * 0.0008
          other_label_y <- y_coords[k] + result$yshift[k] * 0.0008

          # Check for precise geometric overlap
          other_width <- text_widths[k]
          other_height <- text_heights[k]

          other_left <- other_label_x - other_width/2
          other_right <- other_label_x + other_width/2
          other_bottom <- other_label_y - other_height/2
          other_top <- other_label_y + other_height/2

          # Calculate overlap areas
          x_overlap_amount <- max(0, min(label_right, other_right) - max(label_left, other_left))
          y_overlap_amount <- max(0, min(label_top, other_top) - max(label_bottom, other_bottom))

          if(x_overlap_amount > 0 && y_overlap_amount > 0) {
            # Actual overlap - very high penalty
            overlap_area <- x_overlap_amount * y_overlap_amount
            collision_penalty <- collision_penalty + 200 + overlap_area * 1000
          } else {
            # Check proximity
            label_dist <- sqrt((test_x - other_label_x)^2 + (test_y - other_label_y)^2)
            if(label_dist < 0.02) collision_penalty <- collision_penalty + 30
            else if(label_dist < 0.035) collision_penalty <- collision_penalty + 10
          }
        }
      }

      # Calculate total penalty for this distance
      total_penalty <- boundary_penalty + collision_penalty + test_distance * 20

      if(total_penalty < min_penalty) {
        min_penalty <- total_penalty
        best_distance <- test_distance
      }
    }

    # Store the minimum penalty found
    space_scores[j] <- min_penalty

    # Enhanced directional bonuses based on actual available space
    available_space <- .calculate_directional_space(point_x, point_y, angle, plot_margin)
    space_scores[j] <- space_scores[j] - (available_space * 15)  # More generous bonus
  }

  # More intelligent angle selection
  valid_angles <- orthogonal_angles[space_scores < 100]  # Only truly valid positions

  # Sort by quality
  sorted_indices <- order(space_scores)

  # Select best angles more intelligently
  if(length(valid_angles) >= 2) {
    best_angles <- valid_angles[order(space_scores[orthogonal_angles %in% valid_angles])][1:min(2, length(valid_angles))]
  } else if(length(valid_angles) == 1) {
    best_angles <- c(valid_angles, orthogonal_angles[sorted_indices[2]])
  } else {
    # All have high penalties, take the two least bad
    best_angles <- orthogonal_angles[sorted_indices[1:2]]
  }

  return(list(
    best_angles = best_angles,
    all_angles = orthogonal_angles,
    scores = space_scores,
    debug_info = data.frame(angle = angle_names, score = space_scores)
  ))
}

# Calculate actual space available in a specific direction
.calculate_directional_space <- function(point_x, point_y, angle, plot_margin) {
  if(angle == 0) {  # Right
    return(min(1, (1 - plot_margin - point_x) / 0.3))
  } else if(angle == pi/2) {  # Top
    return(min(1, (1 - plot_margin - point_y) / 0.3))
  } else if(angle == pi) {  # Left
    return(min(1, (point_x - plot_margin) / 0.3))
  } else if(angle == 3*pi/2) {  # Bottom
    return(min(1, (point_y - plot_margin) / 0.3))
  }
  return(0)
}

