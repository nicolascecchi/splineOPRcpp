

library(splineOP)
set.seed(1)
K <- 2
segments <- generate_segment_lengths(20, K,alpha = rep(100,K))
accelerations <- c(5,-30)
result <- generate_Qsplines(segments, accelerations, max1 = TRUE)
signal <- generate_Qspline_signal(result, segments, noise_sd = 0.025)
Rcpp::sourceCpp("./src/SplineOP_Rcpp.cpp")

newsignal = rbind(signal,signal)
data_matrix <- matrix(as.numeric(newsignal), nrow=2, ncol=20)

data_matrix
for (seed in 1:25){
  spop <- new(SplineOP, data_matrix
              , 3 # states
              , 3 # speeds
              ,100 # variance
              ,seed # seed
              )
  print(seed)
  print(cumsum(segments))
  spop$predict(1.125)
  print(spop$get_changepoints)
}
spop <- new(SplineOP, data_matrix
            , 3 # states
            , 3 # speeds
            ,100 # variance
            ,seed # seed
)
spop$predict(0.0125)
spop$get_changepoints
spop$get_states
spop$get_argmin_i

#spop$get_speeds[53]
