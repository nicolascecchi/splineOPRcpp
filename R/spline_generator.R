#' Generate Random Segment Lengths Summing to N
#'
#' @description
#' Generates `K` integer segment lengths that sum exactly to `N`, based on
#' proportions drawn from a Dirichlet distribution with parameter vector `alpha`.
#'
#' @param N Integer. The total length (sum of all segment lengths).
#' @param K Integer. The number of segments.
#' @param alpha Numeric vector of length `K`. Parameters of the Dirichlet
#'   distribution. Defaults to a uniform vector of ones.
#'
#' @details
#' The algorithm proceeds as follows:
#' 1. Draws proportions from a Dirichlet distribution using gamma variables.
#' 2. Multiplies by `N` and floors the result to obtain initial segment lengths.
#' 3. Adjusts the lengths so that they sum exactly to `N` by distributing the
#'    remaining difference to the segments with the largest fractional parts.
#'
#' @return
#' An integer vector of length `K` whose elements are non-negative and sum to `N`.
#'
#' @examples
#' # Generate 5 segments summing to 100
#' generate_segment_lengths(100, 5)
#'
#' # Use a non-uniform Dirichlet parameter
#' generate_segment_lengths(50, 3, alpha = c(2, 5, 3))
generate_segment_lengths <- function(N, K, alpha = rep(1, K))
{
  # Internal Dirichlet sampler
  rdirichlet <- function(n, alpha)
  {
    X <- matrix(rgamma(n * length(alpha), alpha), ncol = length(alpha))
    return(as.vector(X / rowSums(X)))
  }

  # Validate inputs
  if (length(alpha) != K)
    stop("Length of 'alpha' must equal K.")
  if (N <= 0 || K <= 0)
    stop("'N' and 'K' must be positive integers.")

  # Draw proportions from a Dirichlet distribution
  proportions <- rdirichlet(1, alpha)
  segment_lengths <- floor(proportions * N)

  # Adjust to ensure the sum equals N
  diff <- N - sum(segment_lengths)
  if (diff != 0)
  {
    # Compute fractional parts (residuals)
    residuals <- (proportions * N) - segment_lengths
    # Order indices by descending residuals
    order_idx <- order(residuals, decreasing = TRUE)
    # Distribute remaining counts to the top |diff| entries
    segment_lengths[order_idx[seq_len(abs(diff))]] <-
      segment_lengths[order_idx[seq_len(abs(diff))]] + sign(diff)
  }
  return(segment_lengths)
}

################################################################################
################################################################################
################################################################################


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
Ksplines <- function(segments, a)
{
  # Input checks
  if (length(segments) != length(a))
  {
    stop("'segments' and 'a' must have the same length.")
  }
  if (any(segments <= 0))
  {
    stop("'segments' must contain only positive values.")
  }

  K <- length(segments)

  # Right-hand side vector
  b <- c(segments * a, 0)

  # Construct coefficient matrix A
  A <- diag(K + 1)
  for (d in seq_len(K) - 1)
  {
    A[row(A) == col(A) + d] <- -2 * (-1)^d
  }

  # Boundary condition (last row)
  A[K + 1, ] <- c(0, segments)

  # Solve the system A * vtilde = b
  vtilde <- solve(A, b)

  # Compute positions and velocities
  positions <- c(0, cumsum(vtilde[-1] * segments))
  speeds <- vtilde[-1] - segments * a

  # Normalize
  scale_factor <- max(abs(positions))
  if (scale_factor > 0)
  {
    positions <- positions / scale_factor
    speeds <- speeds / scale_factor
    a <- a / scale_factor
  }

  return(list(p = positions, v = speeds, a = a))
}


################################################################################
################################################################################
################################################################################


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
#' res <- Ksplines(seg, acc)
#' plot_Ksplines(res$p, res$v, res$a, seg)
#'
#' @export
plot_Ksplines <- function(positions, speeds, accelerations, segment_lengths,
                          main = "Quadratic Spline", col = "blue")
{

  Ksplines_data <- function(positions, speeds, accelerations, segment_lengths, sd = 1)
  {
    y <- NULL
    x_start <- 0
    for (i in 1:(length(seq_along(positions))-1))
    {
      t <- seq(0,segment_lengths[i]-1, length.out = segment_lengths[i])  # Fine sampling
      y <- c(y,positions[i] + speeds[i] * t + accelerations[i] * t^2 + rnorm(n = length(t), sd = sd))
      x_start <- x_start + segment_lengths[i]
    }
    return(y)
  }

  total_length <- sum(segment_lengths)
  plot(NULL, xlim = c(0, total_length), ylim = c(min(y), max(y)),
       xlab = "Time", ylab = "Position", main = main)

  for (i in seq_len(length(positions)))
  {
    t <- seq(0, segment_lengths[i], length.out = 200)
    y <- y_start + speeds[i] * t + 0.5 * accelerations[i] * t^2
    lines(x_start + t, y, col = col, lwd = 2)
    x_start <- x_start + segment_lengths[i]
    if (i < length(positions))
      y_start <- tail(y, 1)  # ensure continuity
  }


  legend("topright", legend = c("Position"), col = col, lty = 1, lwd = 2, bty = "n")
}






