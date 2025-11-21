path_cache_file <- "~/Desktop/path.txt"
if (file.exists(path_cache_file)) {
  Sys.setenv(PATH = readLines(path_cache_file)[1])
} else {
  writeLines(Sys.getenv("PATH"), path_cache_file)
}
Rcpp::sourceCpp("./src/SplineOP_Rcpp.cpp")
library(tidyverse)
library(splineOP)
library(foreach)
library(doParallel)
num_cores <- detectCores() - 1 # Use all but one core for stability
cl <- makeCluster(num_cores)  # Create a cluster object
registerDoParallel(cl)        # Register the cluster for foreach
print(paste("Using", num_cores, "cores for parallel processing."))
clusterEvalQ(cl, {
  library(Rcpp)
  library(splineOP)
})

# 3. Export Custom Functions/Objects
custom_functions <- c("sample_and_flip")
clusterExport(cl, varlist = custom_functions)


# Setup
#for(seed in 1:100)

seeds <- 1:10
nstates_values <- c(5, 10, 15, 20, 25)
N_values <- c(100, 200, 500)
K_values <- c(2,3,5,10,15)

job_list <- expand.grid(
  seed = seeds,
  nstates = nstates_values,
  N = N_values,
  K = K_values
)
num_runs <- nrow(job_list)

library(pbapply)
job_list_rows <- split(job_list, seq(nrow(job_list)))

run_single_job <- function(job_params) {
  Rcpp::sourceCpp("./src/SplineOP_Rcpp.cpp")
  seed <-job_params$seed
  nstates <- job_params$nstates
  N <- job_params$N
  K <- job_params$K

  set.seed(seed)
  coincide <-  0
  segments <- generate_segment_lengths(N, K,alpha = rep(100,K))
  accelerations1 <- sample_and_flip(seed, K)
  result1 <- generate_Qsplines(segments, accelerations1, max1 = TRUE)
  clean1 <- generate_Qspline_signal(result1, segments)
  noised1 <- add_noise(clean1,30)

  # Get noise
  sigma_signal1 <- var(clean1)
  sigma_signal <- sigma_signal1
  noise_sd = 10^(-30/20) * sigma_signal

  # Make it 2D (but repeated)
  data_matrix <- matrix(as.numeric(noised1), nrow=1, ncol=N)

  # data, nstates, nspeeds, speed estimator lenght, variance for state generation, seed for state generation
  sp = as.numeric(c(35))
  spop1 <- new(SplineOP, data_matrix , nstates , 5, sp, noise_sd*sqrt(nstates), seed )
  spop2 <- new(SplineOP, data_matrix , nstates , 5, sp, noise_sd*sqrt(nstates), seed )
  initspeed = matrix(result1$v[[1]])
  spop1$set_initSpeeds(initspeed)
  spop2$set_initSpeeds(initspeed)

  base_beta = 0.00125
  beta = base_beta*2/(K)
  system.time(spop1$pruningv1(beta, 0))
  system.time(spop2$predict(beta)) #3.86sec

  coincide = identical(spop1$get_changepoints,spop2$get_changepoints)

  data.frame(
    seed = seed,
    nstates = nstates,
    N = N,
    K=K,
    coincide = coincide
  )
}

job_list_rows <- split(job_list, seq(nrow(job_list)))
results_list <- pblapply(X = job_list_rows,
                         FUN = run_single_job,
                         cl = cl)

# Combine the results list back into a single data frame
results_parallel <- do.call(rbind, results_list)
#write.csv(results_parallel, "theorem3-20251118-run2.csv")

summary_df <- results_parallel %>%
  group_by(nstates, N, K) %>%
  summarise(
    total_coincide = sum(coincide, na.rm = TRUE),
    .groups = 'drop' # Recommended to drop the grouping structure afterwards
  )

ggplot(summary_df, aes(x = N, y = total_coincide, color = factor(K), shape = factor(nstates))) +
  geom_point(size = 3, alpha = 0.8) + # Use geom_point for a dot plot
  labs(
    title = "Coincidence by N, colored by K, shaped by Nstates",
    x = "N (X-axis)",
    y = "Total Coincidence (Y-axis)",
    color = "K Value",
    shape = "Number of States"
  ) +
  theme_minimal() + # A clean theme for the plot
  theme(legend.position = "right")
# Run the parameter sweep in parallel
results_parallel <- foreach(i = 1:num_runs,
                            .combine = 'rbind',
                            .packages = c('base', 'stats', 'splineOP','Rcpp')) %dopar% {

      Rcpp::sourceCpp("./src/SplineOP_Rcpp.cpp")
      seed <- job_list[i, "seed"]
      nstates <- job_list[i, "nstates"]
      N <- job_list[i, "N"]
      K <- job_list[i, "K"]

      set.seed(seed)
      coincide <-  0
      segments <- generate_segment_lengths(N, K,alpha = rep(100,K))
      accelerations1 <- sample_and_flip(seed, K)
      result1 <- generate_Qsplines(segments, accelerations1, max1 = TRUE)
      clean1 <- generate_Qspline_signal(result1, segments)
      noised1 <- add_noise(clean1,30)

      # Get noise
      sigma_signal1 <- var(clean1)
      sigma_signal <- sigma_signal1
      noise_sd = 10^(-30/20) * sigma_signal

      # Make it 2D (but repeated)
      data_matrix <- matrix(as.numeric(noised1), nrow=1, ncol=N)

      # data, nstates, nspeeds, speed estimator lenght, variance for state generation, seed for state generation
      sp = as.numeric(c(35))
      spop1 <- new(SplineOP, data_matrix , nstates , 5, sp, noise_sd*sqrt(nstates), seed )
      spop2 <- new(SplineOP, data_matrix , nstates , 5, sp, noise_sd*sqrt(nstates), seed )
      initspeed = matrix(result1$v[[1]])
      spop1$set_initSpeeds(initspeed)
      spop2$set_initSpeeds(initspeed)

      beta = 0.00125
      system.time(spop1$pruningv1(beta, 0))

      beta = 0.00125
      system.time(spop2$predict(beta)) #3.86sec

      coincide = identical(spop1$get_changepoints,spop2$get_changepoints)

      data.frame(
        seed = seed,
        nstates = nstates,
        N = N,
        K=K,
        coincide = coincide
      )
    }



# Stop the parallel cluster
stopCluster(cl)

# The final results are in 'results_parallel'
head(results_parallel)

results_df <- data.frame(all_combinations, metric_value = NA)

sample_and_flip <- function(seed, K) {
  set.seed(seed)
  sampled_values <- sample(1:20, size = K, replace = TRUE)
  sign_pattern <- rep(c(1, -1), length.out = K)
  final_vector <- sampled_values * sign_pattern
  return(final_vector)
}

row_index <- 1
coincide = 0
for (T in T_values){
  for (K in K_values){
    for (nstates in nstates_values){
      coincide = 0
      for (seed in seeds){
        set.seed(seed)

        segments <- generate_segment_lengths(T, K,alpha = rep(100,K))
        accelerations1 <- sample_and_flip(seed, K)
        result1 <- generate_Qsplines(segments, accelerations1, max1 = TRUE)
        clean1 <- generate_Qspline_signal(result1, segments)
        noised1 <- add_noise(clean1,30)

        # Get noise
        sigma_signal1 <- var(clean1)
        sigma_signal <- sigma_signal1
        noise_sd = 10^(-30/20) * sigma_signal

        # Make it 2D (but repeated)
        data_matrix <- matrix(as.numeric(noised1), nrow=1, ncol=T)

        # data, nstates, nspeeds, speed estimator lenght, variance for state generation, seed for state generation
        sp = as.numeric(c(35))
        spop1 <- new(SplineOP, data_matrix , nstates , 5, sp, noise_sd*sqrt(nstates), seed )
        spop2 <- new(SplineOP, data_matrix , nstates , 5, sp, noise_sd*sqrt(nstates), seed )
        initspeed = matrix(result1$v[[1]])
        spop1$set_initSpeeds(initspeed)
        spop2$set_initSpeeds(initspeed)

        beta = 0.00125
        system.time(spop1$pruningv1(beta, 0))

        beta = 0.00125
        system.time(spop2$predict(beta)) #3.86sec

        coincide = coincide + identical(spop1$get_changepoints,spop2$get_changepoints)
      }
      print(paste0("K ", K, " nstates ", nstates, " coindide ", coincide))
      results_df[row_index, "metric_value"] <- coincide
      row_index <-row_index+1
    }
  }
}
#write_csv(x = results_df, file = "results_analysis.csv")
results_df



{ seed=2
  print(seed)
  set.seed(seed)
  K <- 5
  T = 200
  segments <- generate_segment_lengths(T, K,alpha = rep(100,K))
  accelerations1 <- sample_and_flip(seed, K)
  result1 <- generate_Qsplines(segments, accelerations1, max1 = TRUE)
  clean1 <- generate_Qspline_signal(result1, segments)
  noised1 <- add_noise(clean1,30)

  # Get noise
  sigma_signal1 <- var(clean1)
  sigma_signal <- sigma_signal1
  noise_sd = 10^(-30/20) * sigma_signal

  # Make it 2D (but repeated)
  data_matrix <- matrix(as.numeric(noised1), nrow=1, ncol=T)

  # data, nstates, nspeeds, speed estimator lenght, variance for state generation, seed for state generation
  sp = as.numeric(c(35))
  NS = 50
  spop1 <- new(SplineOP, data_matrix , NS , 5, sp, noise_sd*NS, 12 )
  spop2 <- new(SplineOP, data_matrix , NS , 5, sp, noise_sd*NS, 12 )
  initspeed = matrix(result1$v[[1]])
  spop1$set_initSpeeds(initspeed)
  spop2$set_initSpeeds(initspeed)

  real = c(1, cumsum(segments))
  beta = 0.00125
  print(system.time(spop1$pruningv1(beta, 0/4)))

  beta = 0.00125
  print(system.time(spop2$predict(beta))) #3.86sec

  coincide = coincide + identical(spop1$get_changepoints,spop2$get_changepoints)
  print(c(1, cumsum(segments)))
  print(spop1$get_changepoints)
  print(spop2$get_changepoints)
}
coincide



# 200 states Precision: Recall:

# Precision and recall metric for change point detection with a margin of error.
#
# This function is adapted from Python's ruptures package, using a margin
# of error to determine a match. It correctly excludes the first and last
# indices (start and end points) from the count of detected changepoints.
#
# @param true_bkps A numeric vector of true indices (includes start and end points).
# @param my_bkps A numeric vector of estimated indices (includes start and end points).
# @param margin An integer specifying the allowed error (in points).
# @return A numeric vector of length 2: (precision, recall).
#
precision_recall <- function(true_bkps, my_bkps, margin = 10) {

  # --- 1. Sanity Checks and Base Case ---

  if (!is.numeric(true_bkps) || !is.numeric(my_bkps)) {
    stop("Inputs 'true_bkps' and 'my_bkps' must be numeric vectors.")
  }

  if (margin < 0) {
    stop(paste0("Margin of error must be positive (margin = ", margin, ")."))
  }

  # If length is <= 2, it means zero actual changepoints were detected (only start and end).
  if (length(my_bkps) <= 2 || length(true_bkps) <= 2) {
    return(c(precision = 0, recall = 0))
  }

  # --- 2. Prepare Vectors (Removing the start and end index) ---

  # Exclude the first index (start) and the last index (end) from both lists.
  true_cp <- true_bkps[-c(1, length(true_bkps))]
  my_cp <- my_bkps[-c(1, length(my_bkps))]

  # --- 3. Core Logic: Vectorized Matching (Replicates Python's product and filter) ---

  # Check 1: Is true_cp within my_cp + margin?
  # Outer computes a matrix where entry [i, j] is my_cp[i] vs true_cp[j]
  upper_bound_check <- outer(my_cp, true_cp, FUN = function(m, t) t < m + margin)

  # Check 2: Is true_cp within my_cp - margin?
  lower_bound_check <- outer(my_cp, true_cp, FUN = function(m, t) t > m - margin)

  # match_matrix is TRUE where an estimated changepoint 'm' covers a true changepoint 't'
  match_matrix <- upper_bound_check & lower_bound_check

  # --- 4. Identify True Positives ---

  # The Python set logic ensures that each *estimated* changepoint (my_cp)
  # is counted as a match only once, even if it falls within the margin of two
  # different true changepoints.
  # In our matrix: rows correspond to estimated points (my_cp).

  # rowSums > 0 means the estimated changepoint at that row found AT LEAST ONE match.
  # This correctly counts the number of estimated changepoints that are "correct".
  matched_my_cp <- rowSums(match_matrix) > 0

  # Count the number of true positives (tp_)
  tp_ <- sum(matched_my_cp)

  # --- 5. Calculate Precision and Recall ---

  num_my_cp <- length(my_cp)
  num_true_cp <- length(true_cp)

  # Calculate metrics, protecting against division by zero (though length check guards this)
  precision <- if (num_my_cp == 0) 0 else tp_ / num_my_cp
  recall <- if (num_true_cp == 0) 0 else tp_ / num_true_cp

  return(c(precision = precision, recall = recall))
}


