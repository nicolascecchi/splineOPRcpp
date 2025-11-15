
path_cache_file <- "~/Desktop/path.txt"
if (file.exists(path_cache_file)) {
  Sys.setenv(PATH = readLines(path_cache_file)[1])
} else {
  writeLines(Sys.getenv("PATH"), path_cache_file)
}
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
