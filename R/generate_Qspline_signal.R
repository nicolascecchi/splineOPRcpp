#' Generate a Quadratic Spline Signal with Additive Noise
#'
#' @param result List output from `generate_Qsplines`, containing vectors `p`, `v`, and `a`
#'   for the position, velocity, and acceleration of each segment.
#' @param segments Numeric vector of segment lengths used to generate the spline.
#' @param n_points Number of points per segment (default 50).
#' @param noise_sd Standard deviation of additive Gaussian noise (default 0).
#'
#' @return A numeric vector representing the continuous spline signal with noise.
#' @examples
#' result <- list(p = c(0, 1, 2), v = c(1, 0.5, -0.5), a = c(0.1, -0.1, 0))
#' segments <- c(3, 2, 4)
#' signal <- generate_Qspline_signal(result, segments, n_points = 100, noise_sd = 0.05)
#' plot(signal, type = "l")
generate_Qspline_signal <- function(result, segments)
{
  x_full <- c()
  p_full <- c()
  cum_len <- c()

  for (i in seq_along(segments))
  {
    seg_len <- segments[i]
    t <- seq(0, seg_len-1, length.out = seg_len)

    # Quadratic interpolation within segment
    p_seg <- result$p[i] + result$v[i] * t + 0.5 * result$a[i] * t^2

    x_full <- c(x_full, cum_len + t)
    p_full <- c(p_full, p_seg)
    cum_len <- cum_len + seg_len
  }

  p_full <- p_full

  return(p_full)
}
