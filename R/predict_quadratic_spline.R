#' Predicts Least-Squares Quadratic Spline
#'
#' Fits a quadratic B-spline (degree 2) to multidimensional data using
#' Ordinary Least Squares (OLS). The spline is fitted to each dimension
#' (column) independently using the same basis.
#'
#' @param data A numeric matrix of shape [N_OBS, N_DIMS].
#' @param knots A numeric vector of internal knots in range [1,NOBS],
#' not including endpoints.
#' @return A numeric matrix of shape [N_OBS, N_DIMS] containing the
#'   fitted values (the smoothed spline curves).
#' @importFrom splines bs
predict_quadratic_spline <- function(data, knots) {
  N <- nrow(data)
  if (max(knots) == 1)
  {
    knots = knots * N
    knots = head(knots,-1)
  }
  # 1. Define the time domain (indices)
  x <- 1:N

  # 2. Generate Quadratic B-Spline Basis
  # degree = 2 (Quadratic)
  # intercept = TRUE (Includes the constant term in the basis)
  B <- splines::bs(x, knots = knots, degree = 2, intercept = TRUE)

  # 3. Solve Least Squares: Y = B * Beta + E
  # We regress the data matrix on the basis.
  # "-1" is used because 'B' already provides the intercept column.
  model <- lm(data ~ B - 1)

  # 4. Return the fitted values (The spline curve)
  return(model$fitted.values)
}
