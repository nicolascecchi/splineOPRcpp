#' Add Gaussian Noise to Signal (1D or ND)
#'
#' Adds independent Gaussian noise to an input signal (vector or matrix) based on
#' a specified Signal-to-Noise Ratio (SNR). Handles both one-dimensional vectors
#' and multi-dimensional matrices/data frames.
#'
#' @param signal A numeric vector (1D) or a matrix/data frame (ND) representing
#'   the input signal. Each dimension is a column, each observation a row.
#' @param snr A numeric value (in dB) specifying the Signal-to-Noise Ratio (SNR).
#'   Defaults to 30 dB.
#' @return The noisy signal, returned as an object of the same shape and type
#'   as the input \code{signal}.
#' @details
#'   - For 1D signals, noise variance is based on the total signal variance.
#'   - For ND signals (matrices), the standard deviation of the signal ($\sigma_{\text{signal}}$)
#'     is calculated from the mean of the variances of the individual dimensions (the diagonal
#'     elements of the covariance matrix).
add_noise <- function(signal,snr=30){
  if (class(signal)=="numeric")
  {
    signal = as.matrix(signal)
  }
  if (ncol(signal)==1)
  {
    sigma_signal <- sqrt(var(signal))
    noise_sd = 10^(-snr/20) * sigma_signal

    noise <- rnorm(length(signal), mean = 0, sd = noise_sd)
    noised_signal = signal + noise
  }else
  {
    varMatrix <- var(signal)
    sigma_signal <- mean(sqrt(diag(varMatrix)))
    noise_sd = 10^(-snr/20) * sigma_signal

    nobs = nrow(signal)
    ndims = ncol(signal)
    noise <- rnorm(nobs*ndims, mean = 0, sd = noise_sd)
    noise_matrix <- matrix(noise, nrow = nobs, ncol = ndims)
    noised_signal = signal + noise_matrix
  }
  return(noised_signal)
}


