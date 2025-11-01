

library(splineOP)
set.seed(51)
K <- 5
segments <- generate_segment_lengths(1000, K,alpha = rep(100,K))
accelerations <- c(5,-3,3,-4,6)
result <- generate_Qsplines(segments, accelerations, max1 = TRUE)
result$p  # positions
result$v  # velocities
result$a  # accelerations
plot_Qspline(result, segments)
signal <- generate_Qspline_signal(result, segments, noise_sd = 0.025)
plot(signal, type = "l")
cumsum(segments)

Rcpp::sourceCpp("./src/SplineOP_Rcpp.cpp")
spop <- new(SplineOP, signal, 5, 1,0.025,5)
spop$predict(0.025)
spop$get_changepoints
cumsum(segments)


spop$get_segment_cost(0,2,0,4,0)
spop$get_segment_cost(2,3,4,3,4.01748)


spop$get_segment_cost(2,4,4,-8,4.01748)
spop$get_segment_cost(2,5,4,-29,4.01748)


spop$predict(250.)
spop$get_costs
spop$get_speeds

spop$get_changepoints

spop$get_argmin_i
spop$get_argmin_s

cost_obj <- new(QuadraticCost, signal)
cost_obj$segmentcost(0,4,0,0,0)
cost_obj$a




set.seed(51)
# Example with 5 segments and given accelerations
K <- 2
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

signal <- generate_Qspline_signal(result, segments, noise_sd = 0.05)
plot(signal, type = "l")
cumsum(segments)

dp_matrix_Rcpp(signal, beta = 0, S = 1, data_var = 0.05)

i = 1
cost_obj$segmentcost(0,i,0,signal[[i]],-0.017485296)


cost_obj$segmentcost(213,400,-1,0,0.008095)
segments

n= 1
L= 1
a= -10
b= 4
c= 4
sum_y= 4
sum_y2= 16
sum_yL1= 0
sum_yL2= 32
sum_L1= 1
sum_L2= 1
sum_L3= 1
sum_L4= 1

cost = 0.0
cost = cost+ a * a * sum_L4 / 4.;
cost = cost + a * b * sum_L3;
cost = cost + (a * c + b * b) * sum_L2;
cost = cost + 2.0 * b * c * sum_L1;
cost = cost + c * c * n;
cost = cost -  a * sum_yL2;
cost = cost - 2.0 * b * sum_yL1;
cost = cost - 2.0 * c * sum_y;
cost = cost + sum_y2;
cost
