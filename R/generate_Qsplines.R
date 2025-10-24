#' Compute Quadratic K-Splines from Segment Lengths and Accelerations
#'
#' @description
#' Computes piecewise quadratic spline quantities (position, speed, acceleration)
#' for a system of segments characterized by given lengths and accelerations.
#'
#' @param segments Numeric vector of positive segment lengths.
#' @param a Numeric vector of accelerations, of the same length as `segments`.
#'
#' @details
#' The function constructs and solves a linear system \eqn{A vtilde = b}, where:
#' \itemize{
#'   \item \eqn{b = (segments * a, 0)} defines acceleration constraints,
#'   \item \eqn{A} is a structured matrix encoding spline continuity conditions,
#'   \item \eqn{v} are auxiliary velocity variables.
#' }
#'
#' The resulting positions, velocities, and accelerations are normalized so that
#' the maximum absolute position equals 1 (unless all positions are zero).
#'
#' @return
#' A list with components:
#' \describe{
#'   \item{p}{Normalized positions along the spline.}
#'   \item{v}{Normalized velocities at segment boundaries.}
#'   \item{a}{Normalized accelerations (same as input, scaled).}
#' }
#'
#' @examples
#' # Example with 5 segments and given accelerations
#' segments <- c(10, 20, 15, 10, 20)
#' accelerations <- c(0.5, -0.2, 0.3, 0.1, -0.1)
#' result <- generate_Qsplines(segments, accelerations)
#'
#' # Inspect results
#' result$p  # positions
#' result$v  # velocities
#' result$a  # normalized accelerations
#'
#' # Plot positions along the spline
#' tim <- rep(1:5, segments)
#' plot(result$p[time], type = "b", main = "Quadratic K-Spline Positions",
#'      xlab = "Segment Index", ylab = "Normalized Position")
#'
generate_Qsplines <- function(segments, a, p0 = 0, pK1 = 0)
{


  K <- length(segments) - 1
  # Construct delta vector for A last row
  delta <- numeric(K + 2)
  delta[1] <- 0.5 * segments[1]              # delta_0
  delta[K+2] <- 0.5 * segments[K+1]          # delta_{K+1}
  if(K >= 1) {
    for(i in 1:K) {
      delta[i+1] <- 0.5 * (segments[i] + segments[i+1])
    }
  }

  ##############################
  # Construct A matrix
  A <- matrix(0, nrow = K+2, ncol = K+2)
  for(i in 1:(K+1)) {
    A[i, i] <- -1
    A[i, i+1] <- 1
  }
  A[K+2, ] <- delta

  # Construct b vector
  b <- c(a * segments, pK1 - p0)

  # Solve the system
  v <- solve(A, b)
  p <- numeric(K+2)
  p[1] <- p0  # initial position
  for(k in 0:(K-1))
  {
    # sum of trapezoids up to segment k
    p[k+2] <- p[k+1] + 0.5 * (v[k+1] + v[k+2]) * delta_tau[k+1]
  }

  return(list(p = p[-length(p)], v = v[-length(v)], a = a))
}












