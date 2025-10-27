

library(splineOP)

# Example with 5 segments and given accelerations
K <- 3
segments <- generate_segment_lengths(400, K,alpha = rep(100,K))
accelerations <- rnorm(K)
result <- generate_Qsplines(segments, accelerations, max1 = TRUE)
#'
# Inspect results
result$p  # positions
result$v  # velocities
result$a  # accelerations
#'
# Plot positions along the spline

plot_Qspline(result, segments)

signal <- generate_Qspline_signal(result, segments, noise_sd = 0.005)
plot(signal, type = "l")


cumsum(segments)
dp_matrix_Rcpp(signal, beta = 200, S = 12, data_var = 0.01)



