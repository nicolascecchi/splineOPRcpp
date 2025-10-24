#' Plot a Quadratic Spline from Segment Parameters
#'
#' @description
#' Plots a piecewise quadratic spline defined by positions, speeds, and accelerations
#' at the start of each segment.
#'
#' @param positions Numeric vector of positions at the start of each segment.
#' @param speeds Numeric vector of initial speeds for each segment.
#' @param accelerations Numeric vector of accelerations for each segment.
#' @param segment_lengths Numeric vector of segment lengths.
#' @param main Character. Title of the plot.
#' @param col Character. Color of the spline curve.
#'
#' @details
#' For each segment \eqn{i}, the spline is defined by:
#' \deqn{y_i(t) = p_i + v_i t + \tfrac{1}{2} a_i t^2, \quad t \in [0, \text{segment\_length}_i].}
#' The function ensures positional continuity across segments.
#'
#' @examples
#' seg <- c(10, 20, 10)
#' acc <- c(0.5, -0.2, 0.3)
#' res <- generate_Qsplines(seg, acc)
#' plot_Qsplines(res$p[-length(res$p)], res$v, res$a, seg)
#'
#' @export
plot_Qsplines <- function(positions,
                          speeds,
                          accelerations,
                          segment_lengths,
                          main = "Quadratic Spline",
                          col = "blue")
{

  if (length(positions) != length(speeds) ||
      length(positions) != length(accelerations) ||
      length(positions) != length(segment_lengths)) {
    stop("All input vectors must have the same length.")
  }

  K <- length(positions)
  x_start <- 0
  y_start <- positions[1]

  # Compute overall min and max for plotting
  y_values <- numeric(0)
  for (i in seq_len(K))
  {
    t <- seq(0, segment_lengths[i], length.out = 200)
    y_seg <- y_start + speeds[i] * t + 0.5 * accelerations[i] * t^2
    y_values <- c(y_values, y_seg)
    y_start <- tail(y_seg, 1)  # ensure continuity
  }

  # Plot setup
  total_length <- sum(segment_lengths)
  plot(NULL, xlim = c(0, total_length), ylim = range(y_values),
       xlab = "Time", ylab = "Position", main = main)

  # Reset for drawing segments
  x_start <- 0
  y_start <- positions[1]

  for (i in seq_len(K))
  {
    t <- seq(0, segment_lengths[i], length.out = 200)
    y_seg <- y_start + speeds[i] * t + 0.5 * accelerations[i] * t^2
    lines(x_start + t, y_seg, col = col, lwd = 2)
    x_start <- x_start + segment_lengths[i]
    y_start <- tail(y_seg, 1)  # update start for next segment
  }

  legend("topright", legend = "Position", col = col, lty = 1, lwd = 2, bty = "n")
}






