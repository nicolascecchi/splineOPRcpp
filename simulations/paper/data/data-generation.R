library(splineOP)

SAVE_DIRECTORY = "/home/cecchi/projects/splineOPRcpp/simulations/paper/data/"

set.seed(12)
N = 1000
for (K_segments in c(5,10, 15)){# 
  for (sample in 1:1000){
    Kfill = sprintf("%0*d", 2, K_segments)
    samplefill = sprintf("%0*d", 2, sample)
    Nfill = sprintf("%0*d", 4, N)

    y <- NULL
    segments <- generate_segment_lengths(N, K_segments, alpha = rep(5, K_segments))
    a <- NULL
    positions <- NULL
    speeds <- NULL
    accelerations <- NULL

    for (i in 1:3){

      # Segment accelerations
      accelerations <- runif(K_segments, min = 10, max = 15)
      accelerations <- accelerations*(-1)**(1:length(accelerations))

      segments <- generate_segment_lengths(N, K_segments,alpha = rep(100,K))
      result <- generate_Qsplines(segments, accelerations, max1 = TRUE)
      sample_dim <- generate_Qspline_signal(result, segments)

      y <- cbind(y, sample_dim)

    }

    # Saving routine
    output_dir = paste0(SAVE_DIRECTORY,"K",K,"/raw/")
    if (!dir.exists(output_dir)) {
      dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)
    }

    write.table(y,
                paste0(output_dir,"K",Kfill,"id", samplefill,"N",Nfill, ".csv"),
                col.names=FALSE,
                row.names=FALSE)
    write.table(cumsum(segments)/N,
                paste0(output_dir,"bkps_K",Kfill,"id", samplefill,"N",Nfill,".csv"),
                col.names=FALSE,
                row.names=FALSE)
  }
}
