path_cache_file <- "~/Desktop/path.txt"
if (file.exists(path_cache_file)) {
  Sys.setenv(PATH = readLines(path_cache_file)[1])
} else {
  writeLines(Sys.getenv("PATH"), path_cache_file)
}
Rcpp::sourceCpp("./src/SplineOP_Rcpp.cpp")

library(splineOP)
# Setup
set.seed(1)
K <- 5
segments <- generate_segment_lengths(10000, K,alpha = rep(100,K))
accelerations1 <- c(5,-30, -15, 12, -17)
result1 <- generate_Qsplines(segments, accelerations1, max1 = TRUE)
clean1 <- generate_Qspline_signal(result1, segments)
noised1 <- add_noise(clean1,30)

accelerations2 <- c(-16,8, -5, 2, -17)
result2 <- generate_Qsplines(segments, accelerations2, max1 = TRUE)
clean2 <- generate_Qspline_signal(result2, segments)
noised2 <- add_noise(clean2,30)

# Get noise
sigma_signal1 <- var(clean1)
sigma_signal2 <- var(clean2)
sigma_signal <- sigma_signal1+sigma_signal2
noise_sd = 10^(-30/20) * sigma_signal

# Make it 2D (but repeated)
data_matrix <- matrix(as.numeric(rbind(noised1,noised2)), nrow=2, ncol=10000)
# Load SplineOP Module

# data, nstates, nspeeds, speed estimator lenght, variance for state generation, seed for state generation
sp = as.numeric(c(200,400,600,800,1000))
#spop1 <- new(SplineOP, data_matrix , 10 , 1 , sp,sigma_signal1,12 )
spop <- new(SplineOP, data_matrix , 5 , 5 , sp,sigma_signal,12 )

beta = 0.00125
system.time(spop$pruningv1(beta, beta/2)) # 10k points, 8sec
spop$get_changepoints # 1  1866  4033  6287  8356 10000
cumsum(segments)      # 1  1839  4070  6288  8331 10000

