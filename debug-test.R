

library(splineOP)
set.seed(51)
K <- 5
segments <- generate_segment_lengths(1000, K,alpha = rep(100,K))
accelerations <- c(5,-3,3,-4,6)
result <- generate_Qsplines(segments, accelerations, max1 = TRUE)
signal <- generate_Qspline_signal(result, segments, noise_sd = 0.025)
Rcpp::sourceCpp("./src/SplineOP_Rcpp.cpp")

newsignal = rbind(signal,signal)
data_matrix <- matrix(as.numeric(newsignal), nrow=2, ncol=1000)
nrow(data_matrix)
ncol(data_matrix)
qc= new(QuadraticCost, data_matrix)
cost = qc$interval_cost(0,6,data_matrix[,0],data_matrix[,6],data_matrix[,2] )
print(cost)
