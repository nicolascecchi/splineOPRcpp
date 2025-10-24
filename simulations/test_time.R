
# Install Rcpp if not already installed
# install.packages("Rcpp")

library(Rcpp)

# Create a large numeric vector
x <- runif(1e8)

# -----------------------------
# 1. Built-in R sum
# -----------------------------
time_builtin <- system.time({
  s_builtin <- sum(x)
})

# -----------------------------
# 2. Handmade R sum (loop)
# -----------------------------
handmade_sum <- function(v) {
  s <- 0
  for (i in seq_along(v)) {
    s <- s + v[i]
  }
  s
}

time_handmade <- system.time({
  s_handmade <- handmade_sum(x)
})

# -----------------------------
# 3. Rcpp sum function
# -----------------------------
cppFunction('
double cpp_sum(NumericVector v) {
  double total = 0;
  int n = v.size();
  for (int i = 0; i < n; ++i) {
    total += v[i];
  }
  return total;
}
')

time_cpp <- system.time({
  s_cpp <- cpp_sum(x)
})


time_cpp2 <- system.time({
  s_cpp2 <- cpp_sum2(x)
})



# -----------------------------
# Compare results
# -----------------------------
cat("Built-in sum:", s_builtin, "\n")
cat("Handmade R sum:", s_handmade, "\n")
cat("Rcpp sum:", s_cpp, "\n\n")
cat("Rcpp sum2:", s_cpp2, "\n\n")

cat("Execution time (seconds):\n")
cat("  Built-in :", time_builtin["elapsed"], "\n")
cat("  Handmade :", time_handmade["elapsed"], "\n")
cat("  Rcpp     :", time_cpp["elapsed"], "\n")
cat("  Rcpp2     :", time_cpp2["elapsed"], "\n")




x <- runif(10)
x
cpp_sum(x)



############


seg <- generate_segment_lengths(N = 100, K = 5)
acc <- rnorm(5, sd = 10)
res <- generate_Qsplines(seg, acc)

res

plot_Qsplines(res$p, res$v, res$a, seg)

############

# Suppose the C++ function is already compiled via sourceCpp
states <- matrix(1:6, nrow = 2, ncol = 3)
data <- c(10, 20, 30)
beta <- 0.5

square_states(data, beta, states)



