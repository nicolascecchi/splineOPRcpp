path_cache_file <- "~/Desktop/path.txt"
if (file.exists(path_cache_file)) {
  Sys.setenv(PATH = readLines(path_cache_file)[1])
} else {
  writeLines(Sys.getenv("PATH"), path_cache_file)
}

library(splineOP)
Rcpp::sourceCpp("./src/SplineOP_Rcpp.cpp")

# Setup
set.seed(1)


K <- 5
N <- 2000
snr <- 30
segments <- generate_segment_lengths(N, K,alpha = rep(100,K))
accelerations1 <- c(5,-5, -5, 0, -7)
result1 <- generate_Qsplines(segments, accelerations1, max1 = TRUE)
clean1 <- generate_Qspline_signal(result1, segments)
noised1 <- add_noise(clean1,snr)

accelerations2 <- c(-10,3, -2, 3, -11)
result2 <- generate_Qsplines(segments, accelerations2, max1 = TRUE)
clean2 <- generate_Qspline_signal(result2, segments)
noised2 <- add_noise(clean2,snr)

# Get noise
sigma_signal1 <- var(clean1)
sigma_signal2 <- var(clean2)
sigma_signal <- sigma_signal1+sigma_signal2
noise_sd = 10^(-snr/20) * sigma_signal

# Make it 2D (but repeated)
data_matrix <- matrix(as.numeric(rbind(noised1,noised2)), nrow=2, ncol=N)
# Load SplineOP Module

# data, nstates, nspeeds, speed estimator lenght, variance for state generation, seed for state generation
sp = as.numeric(c(200,400,600,800,1000))/10
#spop1 <- new(SplineOP, data_matrix , 10 , 1 , sp,sigma_signal1,12 )
spop <-  new(SplineOP,data_matrix, 5, 5, sp, noise_sd/4, 12)
spop2 <- new(SplineOP,data_matrix, 5, 5, sp, noise_sd, 12)

penalty = 0.00125
margin = 1.25
system.time(spop$pruningv1(penalty, margin))
system.time(spop2$predict(penalty)) # 10k points, 8sec
spop$get_changepoints # 1  1866  4033  6287  8356 10000
spop2$get_changepoints
sapply(spop$get_non_pruned_times(), length)
#cumsum(segments)
