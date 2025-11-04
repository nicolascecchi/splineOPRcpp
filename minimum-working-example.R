

library(splineOP)
# Setup
set.seed(1)
K <- 4
segments <- generate_segment_lengths(1000, K,alpha = rep(100,K))
accelerations <- c(5,-30, 10, -2)
result <- generate_Qsplines(segments, accelerations, max1 = TRUE)
clean <- generate_Qspline_signal(result, segments)
noised <- add_noise(clean,30)

# Get noise
sigma_signal <- sqrt(var(clean))
noise_sd = 10^(-30/20) * sigma_signal

# Make it 2D (but repeated)
newsignal = rbind(noised,noised)
data_matrix <- matrix(as.numeric(newsignal), nrow=2, ncol=1000)

# Load SplineOP Module
Rcpp::sourceCpp("./src/SplineOP_Rcpp.cpp")
spop <- new(SplineOP, data_matrix
              , 4 # states
              , 5 # speeds
              ,0.0126 # variance
              ,12 # seed
              )
spop$predict(0.0025)
spop$get_changepoints # 226, 487, 765, 999
cumsum(segments) # 221, 489, 755, 1000
