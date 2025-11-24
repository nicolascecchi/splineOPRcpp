#' Estimate standard deviation of the signal.
#'
#' Estimates the standard deviation of a signal.
#' For multi-dimension case, estimates the sd on each dimension and
#' returns the mean value of those estimations.
#' See: Asymptotically Optimal Difference-Based Estimation of Variance in Nonparametric Regression
#' Peter Hall, J. W. Kay and D. M. Titterington
#' Stable URL https://www.jstor.org/stable/2336990
#' @param data The observed data.
#' @return Estimated standard deviation.
sdHallDiff2 <- function(data)
{
  NOBS = ncol(data)
  NDIM = nrow(data)

  if (NDIM == 1)
  {
    wei <- c(0.1942, 0.2809, 0.3832, -0.8582)
    corrector <- wei[4]^2 +
      (2*wei[4]-wei[3])^2 +
      (-wei[2] + 2*wei[3]-wei[4])^2 +
      (-wei[1] + 2*wei[2]-wei[3])^2 +
      (2*wei[1]-wei[2])^2 +
      wei[1]^2
    m <- length(data)
    z <- - data[3:m] - data[1:(m-2)] +  2*data[2:(m-1)]
    n <- length(z)
    mat <- wei %*% t(z)
    mat[2, -n] <- mat[2, -1]
    mat[3, -c(n-1, n)] <- mat[3, -c(1, 2)]
    mat[4, -c(n-2, n-1, n)] <- mat[4, -c(1, 2, 3)]
    sd <- sqrt(sum(apply(mat[, -c(n-2, n-1, n)], 2, sum)^2)/((n-3)*corrector))
    return(sd)
  }
  # Multi-dimensional case
  # Compute the estimator for each dimension and make the mean
  # We assume the noise is the same on all dimension, we are getting a better
  # estimator.
  else
  {
    sd_cum = 0.0
    for (dim in 1:NDIM)
    {
      sd_cum = sd_cum + sdHallDiff2(data[dim,])
    }
    sd_mean = sd_cum/NDIM
    return(sd_mean)
  }
}

