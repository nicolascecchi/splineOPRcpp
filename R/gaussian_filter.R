#' @title Apply Gaussian Smoothing to a Multi-Dimensional Time Series.
#'
#' @description Generates a 1D Gaussian kernel and applies it as a centered
#' convolution filter to each column (dimension) of a matrix or data frame using
#' the \code{signal::conv()} function. This function ensures the output signal
#' has the same length as the input signal (\code{N_out = N_in}).
#'
#' **Boundary Handling:** To maintain length, the signal is padded at both ends
#' by replicating the first and last values for a length equal to the kernel radius.
#'
#' **Kernel Truncation:** The filter kernel is automatically truncated at
#' \code{+/- 4 * sigma} (hardcoded).
#'
#' @param data A numeric matrix or data frame where rows are time steps and
#' columns are the dimensions/variables (e.g., [100 x 3]).
#' @param sigma The standard deviation of the Gaussian distribution. Higher sigma
#' results in smoother output. Must be > 0.
#' @return A matrix of the exact same dimensions as \code{data} containing the
#' smoothed time series, with no \code{NA} values.
#' @importFrom signal conv
#' @export
gaussian_filter_3d <- function(data, sigma) {

  # 1. Input Validation and Coercion
  if (!is.matrix(data) && !is.data.frame(data)) {
    stop("Input 'data' must be a matrix or a data frame.")
  }
  if (sigma <= 0) {
    stop("Sigma must be positive.")
  }

  ts_matrix <- as.matrix(data)
  n_rows <- nrow(ts_matrix)
  n_cols <- ncol(ts_matrix)

  # 2. Determine Filter Radius (Hardcoded truncation at 4 * sigma)
  truncate <- 4.0
  radius <- ceiling(truncate * sigma)

  if (radius < 1) {
    radius <- 1
  }

  # 3. Generate the Gaussian Kernel Weights

  # The kernel range is from -radius to +radius
  kernel_range <- seq(-radius, radius, by = 1)

  # Calculate Gaussian PDF weights
  kernel_weights <- exp(-(kernel_range^2) / (2 * sigma^2))

  # Normalize the kernel so the weights sum to 1
  kernel_weights <- kernel_weights / sum(kernel_weights)

  # 4. Apply the Convolution (Filtering) to Each Column using signal::conv

  # Initialize an output matrix
  smoothed_matrix <- matrix(NA_real_, nrow = n_rows, ncol = n_cols)

  for (j in 1:n_cols) {
    column_data <- ts_matrix[, j]

    # --- Padding (Replicate Boundary Handling) ---
    # Replicate the first and last value 'radius' times for padding.
    start_val <- column_data[1]
    end_val <- column_data[n_rows]

    pad_start <- rep(start_val, radius)
    pad_end <- rep(end_val, radius)

    column_padded <- c(pad_start, column_data, pad_end)

    # ---------------------------------------------

    # Perform convolution on the padded signal
    # The output length is L_conv = L_padded + L_kernel - 1 = (N+2R) + (2R+1) - 1 = N + 4R
    convolved_column <- signal::conv(column_padded, kernel_weights)

    # Trim and Center:
    # We extract the central N elements, which now correspond to the smoothed
    # original signal, including the padded regions for stable filtering at the edges.
    # The central N elements start at index 2*radius + 1.

    start_idx <- 2 * radius + 1
    end_idx <- start_idx + n_rows - 1

    # Extract the central N points
    smoothed_center <- convolved_column[start_idx:end_idx]

    # Assign the smoothed result to the final matrix
    smoothed_matrix[, j] <- smoothed_center
  }

  # Preserve column names
  colnames(smoothed_matrix) <- colnames(ts_matrix)

  return(smoothed_matrix)
}
