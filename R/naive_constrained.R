#' Constrained Naive Jump Detector
#'
#' Finds the K largest jumps in the signal's acceleration, measured by the
#' L2 norm (magnitude) of the third-order differences.
#'
#' @param signal A numeric matrix or vector representing the input data.
#'   It is expected to be of shape [N_SAMPLES, N_DIMS].
#' @param K An integer specifying the number of changepoints to detect.
#' @return A numeric vector containing the chronologically sorted indices of
#'   the K largest jumps, plus the last index.
#' @details
#'   The function computes the third-order difference of the signal (along the
#'   sample axis), calculates the magnitude (L2 norm) of the resulting vector
#'   for each time point, and identifies the indices corresponding to the
#'   K largest magnitudes.
naive_constrained <- function(signal, K)
{
  N_SAMPLES = nrow(signal)
  differences <- diff(signal, differences = 3)
  deltas <- sqrt(rowSums(differences^2))
  sorted_idxs <- order(deltas, decreasing = TRUE)[1:K]
  max_idxs <- sort(sorted_idxs)
  max_idxs <- c(max_idxs, N_SAMPLES)
  return(max_idxs)
}

