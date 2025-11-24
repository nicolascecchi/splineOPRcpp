#' Sample Random Integers and Flip Signs
#'
#' Samples K random integers from the range 1 to 20 (inclusive) and then
#' applies an alternating sign pattern (+1, -1, +1, -1, ...) to the resulting vector,
#' ensuring the first element is positive. The process is made reproducible
#' using the provided seed.
#'
#' @param seed An integer used to initialize the random number generator, ensuring
#'   reproducibility of the sampled values.
#' @param K An integer specifying the number of samples to generate (the length
#'   of the resulting vector).
#' @return A numeric vector of length K containing sampled integers (1 to 20)
#'   with alternating signs, starting with a positive value.
#' @examples
#' # Example 1: Basic usage
#' sample_and_flip(seed = 10, K = 5)
#'
#' # Example 2: Reproducibility check
#' v1 <- sample_and_flip(seed = 42, K = 8)
#' v2 <- sample_and_flip(seed = 42, K = 8)
#' all.equal(v1, v2) # Should be TRUE
sample_alternating_accelerations <- function(K, seed=NULL) {
  if (is.null(seed))
  {
    sampled_values <- sample(1:20, size = K, replace = TRUE)
    sign_pattern <- rep(c(1, -1), length.out = K)
    final_vector <- sampled_values * sign_pattern
    return(final_vector)
  }
  else
  {
    set.seed(seed)
    sampled_values <- sample(1:20, size = K, replace = TRUE)
    sign_pattern <- rep(c(1, -1), length.out = K)
    final_vector <- sampled_values * sign_pattern
    return(final_vector)
  }
}
