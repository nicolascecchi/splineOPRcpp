#' Compute Quadratic Splines from Segment Lengths and Accelerations
#'
#' @description
#' Computes piecewise quadratic spline quantities (positions, velocities, accelerations)
#' for a system of segments characterized by given lengths and accelerations.
#' The spline starts and ends at specified positions (`p0` and `pK1`), which default to zero.
#' Optionally, the positions, velocities, and accelerations can be normalized so that
#' the maximum absolute position equals 1.
#'
#' @param segments Numeric vector of positive segment lengths.
#' @param a Numeric vector of accelerations, of the same length as `segments`.
#' @param max1 Logical, if TRUE the positions, velocities, and accelerations are scaled
#'        so that the maximum absolute position equals 1 (default TRUE).
#' @param p0 Numeric, starting position (default 0).
#' @param pK1 Numeric, ending position (default 0).
#'
#' @details
#' The function constructs and solves a linear system \eqn{A v = b}, where:
#' \itemize{
#'   \item \eqn{b = (segments * a, p_{end} - p_{start})} encodes the acceleration and boundary constraints,
#'   \item \eqn{A} is a structured matrix encoding spline continuity conditions,
#'   \item \eqn{v} are auxiliary velocity variables at segment boundaries.
#' }
#'
#' @return A list with components:
#' \describe{
#'   \item{p}{Positions along the spline (normalized if `max1 = TRUE`).}
#'   \item{v}{Velocities at segment boundaries (normalized if `max1 = TRUE`).}
#'   \item{a}{Accelerations (scaled if `max1 = TRUE`).}
#' }
#'
#' @examples
#' # Example with 5 segments and given accelerations
#' segments <- generate_segment_lengths(400, 5)
#' accelerations <- rnorm(5)
#' result <- generate_Qsplines(segments, accelerations, max1 = TRUE)
#'
#' # Inspect results
#' result$p  # positions
#' result$v  # velocities
#' result$a  # accelerations
#'
#' # Plot positions along the spline
#' plot_Qspline(result, segments)
generate_Qsplines <- function(segments, a, max1 = TRUE, p0 = 0, pK1 = 0)
{
  K <- length(segments) - 1

  # Construct delta vector for last row of A
  delta <- numeric(K + 2)
  delta[1] <- 0.5 * segments[1]
  delta[K + 2] <- 0.5 * segments[K + 1]
  if (K >= 1)
  {
    for (i in 1:K)
    {
      delta[i + 1] <- 0.5 * (segments[i] + segments[i + 1])
    }
  }

  # Construct A matrix
  A <- matrix(0, nrow = K + 2, ncol = K + 2)
  for (i in 1:(K + 1))
  {
    A[i, i] <- -1
    A[i, i + 1] <- 1
  }
  A[K + 2, ] <- delta

  # Construct b vector
  b <- c(a * segments, pK1 - p0)

  # Solve the system
  v <- solve(A, b)

  # Compute positions
  p <- numeric(K + 2)
  p[1] <- p0
  for (k in 0:(K - 1))
  {
    p[k + 2] <- p[k + 1] + 0.5 * (v[k + 1] + v[k + 2]) * segments[k + 1]
  }


  # Exclude the last boundary value
  p_out <- p[-length(p)]
  v_out <- v[-length(v)]
  a_out <- a

  # Optional normalization
  if ((max1 == TRUE) && (max(abs(p_out)) > 0))
  {
    scale_factor <- max(abs(p_out))
    p_out <- p_out / scale_factor
    v_out <- v_out / scale_factor
    a_out <- a_out / scale_factor
  }

  list(
    p = p_out,
    v = v_out,
    a = a_out
  )
}
