#' Calculate Precision and Recall for Changepoint Detection
#'
#' This function computes the precision and recall of an estimated segmentation
#' compared to the true segmentation, allowing for a margin of error. It
#' automatically excludes the first and last elements (start and end points)
#' from the changepoint count.
#'
#' @param true_bkps A numeric vector of true indices (includes start and end points).
#' @param my_bkps A numeric vector of estimated indices (includes start and end points).
#' @param margin An integer specifying the allowed error (in points). Defaults to 10.
#' @return A numeric vector of length 2 with named elements: \code{precision} and \code{recall}.
#' @export
#' @examples
#' # Example: True changepoints at 50 and 100, plus start/end
#' true_cp <- c(1, 50, 100, 200)
#' # Estimated changepoints are slightly off, plus an extra point
#' my_cp <- c(1, 48, 105, 150, 200)
#' # Calculate with a margin of 10
#' precision_recall(true_cp, my_cp, margin = 10)
precision_recall <- function(true_bkps, my_bkps, margin = 10) {

  if (!is.numeric(true_bkps) || !is.numeric(my_bkps)) {
    stop("Inputs 'true_bkps' and 'my_bkps' must be numeric vectors.")
  }

  if (margin <= 0) {
    stop(paste0("Margin of error must be positive (margin = ", margin, ")."))
  }

  if (length(my_bkps) <= 2 || length(true_bkps) <= 2) {
    return(c(precision = 0, recall = 0))
  }

  true_cp <- true_bkps[-c(1, length(true_bkps))]
  my_cp <- my_bkps[-c(1, length(my_bkps))]

  upper_bound_check <- outer(my_cp, true_cp, FUN = function(m, t) t < m + margin)

  lower_bound_check <- outer(my_cp, true_cp, FUN = function(m, t) t > m - margin)

  match_matrix <- upper_bound_check & lower_bound_check

  matched_my_cp <- rowSums(match_matrix) > 0

  tp_ <- sum(matched_my_cp)

  num_my_cp <- length(my_cp)
  num_true_cp <- length(true_cp)

  precision <- if (num_my_cp == 0) 0 else tp_ / num_my_cp
  recall <- if (num_true_cp == 0) 0 else tp_ / num_true_cp

  return(c(precision = precision, recall = recall))
}
