

library(splineOP)
set.seed(5)
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

#signal <- generate_Qspline_signal(result, segments, noise_sd = 0.05)
plot(signal, type = "l")
cumsum(segments)



dp_matrix_Rcpp(signal, beta = 0, S = 1, data_var = 0.05)

signal <- c(0,1,4,3,-8,-29)
Rcpp::sourceCpp("./src/QuadraticCost_Rcpp.cpp")
cost_obj <- new(QuadraticCost, signal)
cost_obj$segmentcost(0,2,0,4,4
                     )


n= 2
L= 2
a= 2
b= 0
c= 0
sum_y= 1
sum_y2= 1
sum_yL1= 1
sum_yL2= 1
sum_L1 =1
sum_L2 =1
sum_L3 =1
sum_L4 =1

# 1. Initialize the variable
cost <- 0.0
# 2. Accumulate the values
cost <- cost + sum_y2
cost <- cost + (a * a * sum_L4)
cost <- cost - (2.0 * a * sum_yL2)
cost
