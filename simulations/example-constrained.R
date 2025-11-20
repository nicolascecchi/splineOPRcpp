path_cache_file <- "~/Desktop/path.txt"
if (file.exists(path_cache_file)) {
  Sys.setenv(PATH = readLines(path_cache_file)[1])
} else {
  writeLines(Sys.getenv("PATH"), path_cache_file)
}
sample_and_flip <- function(K, seed=NULL) {
  if (!is.null(seed)) {
    set.seed(seed)
  }
  sampled_values <- sample(1:20, size = K, replace = TRUE)
  sign_pattern <- rep(c(1, -1), length.out = K)
  final_vector <- sampled_values * sign_pattern
  return(final_vector)
}

Rcpp::sourceCpp("./src/SplineOP_Rcpp.cpp")
library(splineOP)
# Setup
seed = 5
set.seed(seed)
Kseg <- 6
N <- 1000
segments <- generate_segment_lengths(N, Kseg,alpha = rep(100,Kseg))
accelerations1 <- sample_and_flip(K=Kseg, seed=seed)
result1 <- generate_Qsplines(segments, accelerations1, max1 = TRUE)
clean1 <- generate_Qspline_signal(result1, segments)
noised1 <- add_noise(clean1,30)
# Get noise
sigma_signal1 <- var(clean1)
sigma_signal <- sigma_signal1
noise_sd = 10^(-30/20) * sigma_signal

# Make it 2D (but repeated)
newsignal = rbind(noised1)
data_matrix <- matrix(as.numeric(newsignal), nrow=1, ncol=N)
sp = as.numeric(c(25,50,75,100,125))
# Load SplineOP Module
spop <- new(SplineOP_constrained, data_matrix , 5, sp,noise_sd,Kseg-1,seed)
system.time(spop$predict(Kseg-1))
spop$get_changepoints # 153  311  479  677  831 1000
cumsum(segments) # 153  310  478  674  831 1000
