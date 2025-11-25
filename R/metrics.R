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

#' Calculate Mean Squared Error (MSE)
#'
#' Computes the Mean Squared Error between the observed signal and the predicted
#' signal (e.g., the segmented fit).
#'
#' @param observed A numeric vector representing the true signal values.
#' @param predicted A numeric vector representing the predicted signal values (fit).
#' @return A single numeric value: the Mean Squared Error (MSE).
#' @export
#' @examples
#' # Simulate a true signal and a noisy prediction
#' true_signal <- 1:10
#' pred_signal <- c(1.1, 1.9, 3.2, 4.1, 5, 6, 7, 8.5, 9, 10.1)
#' mse(true_signal, pred_signal)
mse <- function(observed, predicted) {

  if (length(observed) != length(predicted)) {
    stop(paste0("Input vectors 'observed' and 'predicted' must have the same length. ",
                "Observed length: ", length(observed),
                ", Predicted length: ", length(predicted), "."))
  }

  # MSE = 1/N * sum((Observed - Predicted)^2)
  squared_errors <- (observed - predicted)^2

  return(mean(squared_errors))
}

#' Calculate F1 Score (Harmonic Mean of Precision and Recall)
#'
#' Computes the F1 Score, which is the harmonic mean of precision and recall.
#' The F1 score is a balanced metric useful for assessing classification performance.
#'
#' @param precision The precision value (e.g., from the \code{precision_recall} function).
#' @param recall The recall value (e.g., from the \code{precision_recall} function).
#' @return A single numeric value: the F1 Score.
#' @export
#' @examples
#' # Example F1 score calculation
#' current_precision <- 0.85
#' current_recall <- 0.70
#' f_score(current_precision, current_recall)
f_score <- function(precision, recall) {

  # Check for valid inputs
  if (is.null(precision) || is.null(recall) || !is.numeric(precision) || !is.numeric(recall)) {
    stop("Precision and recall must be single numeric values.")
  }

  # F1 Score is 2 * (Precision * Recall) / (Precision + Recall)

  if ((precision + recall) == 0) {
    # Handles the edge case where both are 0, resulting in an F1 of 0.
    return(0)
  }

  f1 <- 2 * (precision * recall) / (precision + recall)

  return(f1)
}
