

library(splineOP)
# Setup
set.seed(1)
K <- 4
segments <- generate_segment_lengths(1000, K,alpha = rep(100,K))
accelerations1 <- c(5,-30, 10, -2)
result1 <- generate_Qsplines(segments, accelerations1, max1 = TRUE)
clean1 <- generate_Qspline_signal(result1, segments)
noised1 <- add_noise(clean1,30)

accelerations2 <- c(-1,5,-1, 2)
result2 <- generate_Qsplines(segments, accelerations2, max1 = TRUE)
clean2 <- generate_Qspline_signal(result2, segments)
noised2 <- add_noise(clean2,30)


# Get noise
sigma_signal1 <- var(clean1)
sigma_signal2 <- var(clean2)
sigma_signal <- sqrt(sigma_signal1 + sigma_signal2)
noise_sd = 10^(-30/20) * sigma_signal

# Make it 2D (but repeated)
newsignal = rbind(noised1,noised2)
data_matrix <- matrix(as.numeric(newsignal), nrow=2, ncol=1000)

# Load SplineOP Module
Rcpp::sourceCpp("./src/SplineOP_Rcpp.cpp")
spop <- new(SplineOP, data_matrix
              , 4 # states
              , 5 # speeds
              ,sigma_signal # variance
              ,12 # seed
              )

spop$predict(0.0025)
spop$get_changepoints # 217, 491, 765, 999
cumsum(segments) # 221, 489, 755, 1000
