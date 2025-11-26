# --- 1. Helper Function: Create 3rd Difference Matrix D ---

#' @title Create the 3rd Order Difference Matrix for Trend Filtering
#' @description Generates the sparse matrix D used in 3rd-order L1 trend filtering.
#'   This matrix is an approximation of the third-order discrete derivative operator.
#' @param n Integer, the number of time points (rows in the original data 'y').
#' @return A sparse matrix of class \code{dgCMatrix} (from the \code{Matrix} package)
#'   with dimensions \eqn{(n-3) \times n}.
#' @importFrom Matrix sparseMatrix
#' @export
create_3rd_diff_matrix <- function(n) {
  if (n < 4) {
    stop("Data must have at least 4 time points (n >= 4) for 3rd-order trend filtering.")
  }

  # Indices for the sparse matrix creation
  num_rows <- n - 3
  i <- 1:num_rows

  # The non-zero diagonals are: [1, -3, 3, -1]
  # Corresponding column offsets (j indices): [0, 1, 2, 3] relative to i
  j1 <- i
  j2 <- i + 1
  j3 <- i + 2
  j4 <- i + 3

  # Use the Matrix package to construct the sparse difference matrix D.
  D <- Matrix::sparseMatrix(
    # Row indices (i)
    i = c(i, i, i, i),
    # Column indices (j)
    j = c(j1, j2, j3, j4),
    # Values (x) corresponding to the diagonals
    x = c(rep(1, num_rows), rep(-3, num_rows), rep(3, num_rows), rep(-1, num_rows)),
    # Dimensions
    dims = c(num_rows, n)
  )
  return(D)
}

# --- 2. Core Function: Solve L1 Trend Filtering Problem ---

#' @title Solve Multidimensional L1 Trend Filtering
#' @description Solves the 3rd-order L1 trend filtering problem for multivariate time series data
#'   using convex optimization via the \code{CVXR} package. The objective is to
#'   minimize \eqn{0.5 \cdot ||Y - X||_F^2 + \lambda \cdot ||D \cdot X||_1}, where \eqn{D}
#'   is the third-order difference operator and \eqn{\lambda} is the penalty.
#'
#' @param y Data matrix of shape [n x d] (time points x dimensions).
#' @param penalty Regularization parameter (lambda) for the L1 norm.
#' @return A list containing the optimization result (X, the filtered signal),
#'         the CVXR problem object, and the execution time.
#'         \itemize{
#'           \item \code{tf_prediction}: The filtered signal matrix \eqn{X} of shape \eqn{n \times d}.
#'           \item \code{prob}: The solved \code{CVXR} problem object.
#'           \item \code{runtime}: The time taken to solve the problem, in seconds.
#'           \item \code{status}: The solver status (e.g., "optimal").
#'         }
#' @importFrom CVXR Variable sum_squares p_norm Minimize Problem solve getValue
#' @export
multidim_l1_tf_fit <- function(y, penalty) {
  # Check inputs
  if (!is.matrix(y)) {
    stop("Input 'y' must be a matrix.")
  }

  # Dimensions of the data
  n <- nrow(y) # number of time points (rows)
  d <- ncol(y) # number of dimensions (columns)

  # Create the 3rd Difference Matrix D
  D <- create_3rd_diff_matrix(n)

  # Define the CVXR variable X (the filtered signal)
  X <- CVXR::Variable(n, d)

  # Formulate the objective function:
  # min_X 0.5 * ||Y - X||_F^2 + penalty * ||D %*% X||_1
  # 1. Fidelity term: 0.5 * sum_squares(Y - X)
  fidelity_term <- 0.5 * CVXR::sum_squares(y - X)

  # 2. Regularization term: penalty * norm(D @ X, 1)
  regularization_term <- penalty * CVXR::p_norm(D %*% X, 1)

  # Combine objective and define the problem
  # Note the use of CVXR::Minimize and CVXR::Problem
  objective <- CVXR::Minimize(fidelity_term + regularization_term)
  problem <- CVXR::Problem(objective)

  # Note the use of CVXR::solve
  solution <- CVXR::solve(problem, solver = "ECOS", verbose = FALSE)

  # Check solver status
  if (solution$status != "optimal") {
    warning(paste("Solver did not converge! Status:", solution$status))
    X_prediction <- NA
  } else {
    # Note the use of CVXR::getValue
    y_hat <- solution$getValue(X)
  }

  # Return results
  return(y_hat)
}
