#' Plot a Quadratic Spline from Generated Positions, Velocities, and Accelerations
#'
#' @param result List output from `generate_Qsplines`, containing p, v, a.
#' @param segments Numeric vector of segment lengths used to generate the spline.
#' @param n_points Number of points per segment for smooth plotting (default 50).
#' @param main Title of the plot (default "Quadratic K-Spline Positions").
#'
#' @return A plot of the quadratic spline and segment boundary positions.
plot_Qspline <- function(result,
                         segments,
                         n_points = 50,
                         main = "Quadratic K-Spline Positions")
{
  x_full <- c()  # cumulative distance along the spline
  p_full <- c()  # interpolated positions

  cum_len <- 0
  for (i in seq_along(segments))
  {
    seg_len <- segments[i]
    t <- seq(0, seg_len, length.out = n_points)
    # Quadratic interpolation within segment
    p_seg <- result$p[i] + result$v[i] * t + 0.5 * result$a[i] * t^2
    x_full <- c(x_full, cum_len + t)
    p_full <- c(p_full, p_seg)
    cum_len <- cum_len + seg_len
  }

  # Plot the spline
  plot(x_full, p_full, type = "l", lwd = 2,
       main = main,
       xlab = "Distance along spline",
       ylab = "Normalized Position")

  # Mark segment boundaries
  seg_pos <- cumsum(c(0, segments[-length(segments)]))
  points(seg_pos, result$p, col = "red", pch = 19)
  #legend("topleft", legend = c("Quadratic spline", "Segment boundary positions"),
  #       col = c("black", "red"), lty = c(1, NA), pch = c(NA, 19))
}
