

library(splineOP)
set.seed(51)
K <- 5
segments <- generate_segment_lengths(1000, K,alpha = rep(100,K))
accelerations <- c(5,-3,3,-4,6)
result <- generate_Qsplines(segments, accelerations, max1 = TRUE)
signal <- generate_Qspline_signal(result, segments, noise_sd = 0.025)
Rcpp::sourceCpp("./src/SplineOP_Rcpp.cpp")

newsignal = signal #rbind(signal,signal)
data_matrix <- matrix(as.numeric(newsignal), nrow=1, ncol=1000)


spop <- new(SplineOP, data_matrix
            , 3 # states
            , 3 # speeds
            ,0.025 # variance
            ,5 # seed
            )
cumsum(segments)
spop$predict(2.0)
spop$get_changepoints
spop$get_initSpeeds
result$v
