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
#'
#' # Almost uniform segment lengths
#' generate_segment_lengths(200, 4, alpha = rep(10^3,4))
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
