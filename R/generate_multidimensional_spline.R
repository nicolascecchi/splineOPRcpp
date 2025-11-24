#' Generate a Multidimensional Quadratic Spline Signal with Additive Noise
#'
#' @param ndims Number of dimensions.
#' @param n_points Signal length.
#' @param snr Signal to noise ratio.
#' @param seed Random seed.
#'
#' @return A numeric vector representing the continuous spline signal with noise.
#' @examples
#' result <- list(p = c(0, 1, 2), v = c(1, 0.5, -0.5), a = c(0.1, -0.1, 0))
#' segments <- c(3, 2, 4)
#' signal <- generate_Qspline_signal(result, segments, n_points = 100, noise_sd = 0.05)
#' plot(signal, type = "l")
generate_multidimensional_spline <- function(ndims,Kseg,n_points,snr,seed=NULL)
{
  if (seed != NULL){set.seed(seed)}
  segments <- generate_segment_lengths(n_points, Kseg,alpha = rep(100,Kseg))

  for(dim in 1:ndims)
  {
    accelerations <- sample_and_flip(K=Kseg, seed=seed)
    result <- generate_Qsplines(segments, accelerations, max1 = TRUE)
    clean <- generate_Qspline_signal(result, segments)
    noised <- add_noise(clean,snr)

  }


}
