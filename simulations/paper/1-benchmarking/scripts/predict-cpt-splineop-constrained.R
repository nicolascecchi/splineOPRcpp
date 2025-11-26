# Script for predicting fixed changepoints with SplineOP on simulated signals
# Predicts one signal at a time,
# saves the resulting metrics in a file in ../results
#


#!/usr/bin/env Rscript
library(optparse) # For parsing arguments on command line
library(splineOP) # Our package


# 1. Define Options
option_list = list(
  make_option(c("--N_SAMPLES"), type="integer", default=1000, help="Number of samples"),
  make_option(c("--N_SEG"), type="integer", default=5, help="K value"),
  make_option(c("--snr"), type="double", default=30, help="Signal to Noise Ratio"),
  make_option(c("--sample"), type="integer", default=6, help="Sample ID"),
  make_option(c("--nstates"), type="integer", default=5, help="Number of states"),
  make_option(c("--N_CHANGEPOINTS_PRED"), type="integer", default=25, help="Max K"),
  make_option(c("--OUT_FOLDER"), type="character", default="output", help="Output directory")
)

# 2. Parse Arguments
opt_parser = OptionParser(option_list=option_list)
opt = parse_args(opt_parser)

# 3. Define Main Function
main <- function(args) {
  algorithm <- "splineop-constrained"
  # Set up for saving the result
  cmd_args <- commandArgs(trailingOnly = FALSE)
  needle <- "--file="
  # 1. Find the argument that contains the script file name
  match_arg <- grep(needle, cmd_args)
  # Extract the path from the --file= argument
  CURR_FOLDER = dirname(normalizePath(sub(needle, "", cmd_args[match_arg])))

  #CURR_FOLDER = getwd()
  EXPERIENCE_FOLDER = dirname(CURR_FOLDER) # eg. 1-benchmarking
  SAVE_FOLDER = paste0(EXPERIENCE_FOLDER, "/results/",algorithm,"/", args$OUT_FOLDER)
  if (!dir.exists(SAVE_FOLDER)) {
    dir.create(SAVE_FOLDER, recursive = TRUE)
  }
  # setup variables for file name
  Kfill <- sprintf("%0*d", 2, args$N_SEG)
  samplefill <- sprintf("%0*d", 2, args$sample)
  Nfill <- sprintf("%0*d", 4, args$N_SAMPLES)

  # To load data
  DATA_FOLDER <- paste0(dirname(EXPERIENCE_FOLDER),"/data/synth")
  noised_signal <-  read.csv(paste0(DATA_FOLDER,"/K",args$N_SEG,"/noised/","K",Kfill,"id", samplefill,"N",Nfill,"SNR",args$snr, ".csv"),header=FALSE)
  # Set up parameters for creating SplineOP object
  noised_signal <- as.matrix(sapply(noised_signal,as.numeric))
  noised_signal <- t(noised_signal) # transpose to have the required format for spop
  estimated_std_dev <- sdHallDiff2(noised_signal)
  len_speed_estimators <- as.integer(c(20, 40, 60, 80, 100))
  n_changepoints_pred <- as.integer(args$N_CHANGEPOINTS_PRED)
  states_seed <- as.integer(args$sample)
  nb_of_states <- as.integer(args$nstates)

  # Predict pipeline
  spop <- new(SplineOP_constrained
              ,noised_signal # Data to fit
              ,nb_of_states
              ,len_speed_estimators # Vector fo speeds estimation
              ,estimated_std_dev # For state generation
              ,n_changepoints_pred # nb of changepoints, one less than segments
              ,states_seed # for state generation
              )
  
  time_data <- system.time(spop$predict(n_changepoints_pred))

  # METRICS

  # Changepoint retrieval metrics 
  real_changepoints <- read.csv(paste0(DATA_FOLDER,"/K",args$N_SEG,"/raw/","/bkps_K",Kfill,"id", samplefill,"N",Nfill, ".csv"),header=FALSE)
  real_changepoints <- unlist(real_changepoints$V1 * args$N_SAMPLES)
  real_changepoints <- c(1, real_changepoints)
  pred_changepoints <- spop$get_changepoints
  pred_changepoints <- c(1, pred_changepoints)
  
  precrecall1pct <- precision_recall(real_changepoints, pred_changepoints, margin=0.01*args$N_SAMPLES)
  precrecall2pct <- precision_recall(real_changepoints, pred_changepoints, margin=0.02*args$N_SAMPLES)
  precrecall25pct <- precision_recall(real_changepoints, pred_changepoints, margin=0.025*args$N_SAMPLES)
  precrecall5pct <- precision_recall(real_changepoints, pred_changepoints, margin=0.05*args$N_SAMPLES)

  fscore1pct <- f_score(precrecall1pct['precision'], precrecall1pct['recall'])
  fscore2pct <- f_score(precrecall2pct['precision'], precrecall2pct['recall'])
  fscore25pct <- f_score(precrecall25pct['precision'], precrecall25pct['recall'])
  fscore5pct <- f_score(precrecall5pct['precision'], precrecall5pct['recall'])

  # Signal reconstruction (MSE)
  clean_signal <- read.csv(paste0(DATA_FOLDER,"/K",args$N_SEG,"/raw/","/K",Kfill,"id", samplefill,"N",Nfill, ".csv"),header=FALSE)
  clean_signal <- as.matrix(clean_signal)
  # "untranspose" the noised_signal to long format
  predicted_signal <- predict_quadratic_spline(t(noised_signal), pred_changepoints[-c(1, length(pred_changepoints))])
  mse <- mse(clean_signal, predicted_signal)
  # Save the results 
  FILE_NAME = paste0("spop_constrained_predictions_N",Nfill,
                     "-K",Kfill,
                     "-SNR",args$snr,
                     "-ID",samplefill,
                     "-KPRED",args$N_CHANGEPOINTS_PRED,
                     ".json")
  #SAVE_PATH = paste0(SAVE_FOLDER, FILE_NAME)
  OUT_FILE_PATH <- paste0(SAVE_FOLDER, "/", FILE_NAME)
  
  results_list <- as.list(args)
  # Add the main result vector
  results_list$changepoints <- pred_changepoints
  results_list$real_changepoints <- real_changepoints
  results_list$len_speed_estimators <- len_speed_estimators
  results_list$estimated_std_dev <- estimated_std_dev
  results_list$n_changepoints_pred <- n_changepoints_pred
  results_list$states_seed <- states_seed
  results_list$execution_time_seconds <- time_data["elapsed"]
  results_list$cpu_user_time_seconds <- time_data["user.self"]
  results_list$precision01 <- precrecall1pct["precision"]
  results_list$precision02 <- precrecall2pct["precision"]
  results_list$precision025 <- precrecall25pct["precision"]
  results_list$precision05 <- precrecall5pct["precision"]
  results_list$recall01 <- precrecall1pct["recall"]
  results_list$recall02 <- precrecall2pct["recall"]
  results_list$recall025 <- precrecall25pct["recall"]
  results_list$recall05 <- precrecall5pct["recall"]
  results_list$fscore1pct <- fscore1pct
  results_list$fscore2pct <- fscore2pct
  results_list$fscore25pct <- fscore25pct
  results_list$fscore5pct <- fscore5pct
  results_list$mse <- mse
  results_list$algorithm <- algorithm
  jsonlite::write_json(results_list, OUT_FILE_PATH, pretty = TRUE, unbox=TRUE)
}

# 4. Execute
main(opt)
# system("")
#  Rscript simulations/paper/1-benchmarking/scripts/predict-cpt-splineop-constrained.R --N_SAMPLES 1000 --N_SEG 5 --snr 30 --sample 1 --nstates 5 --N_CHANGEPOINTS_PRED 5 --OUT_FOLDER '20251124'
