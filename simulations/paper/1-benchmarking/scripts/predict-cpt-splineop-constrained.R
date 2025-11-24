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
  make_option(c("--K"), type="integer", default=5, help="K value"),
  make_option(c("--snr"), type="double", default=30, help="Signal to Noise Ratio"),
  make_option(c("--sample"), type="integer", default=6, help="Sample ID"),
  make_option(c("--nstates"), type="integer", default=5, help="Number of states"),
  make_option(c("--KMAX"), type="integer", default=25, help="Max K"),
  make_option(c("--OUT_FOLDER"), type="character", default="output", help="Output directory")
)

# 2. Parse Arguments
opt_parser = OptionParser(option_list=option_list)
opt = parse_args(opt_parser)

# 3. Define Main Function
main <- function(args) {
  # Set up for saving the result
  cmd_args <- commandArgs(trailingOnly = FALSE)
  needle <- "--file="
  # 1. Find the argument that contains the script file name
  match_arg <- grep(needle, cmd_args)
  # Extract the path from the --file= argument
  CURR_FOLDER = dirname(normalizePath(sub(needle, "", cmd_args[match_arg])))

  #CURR_FOLDER = getwd()
  EXPERIENCE_FOLDER = dirname(CURR_FOLDER) # eg. 1-benchmarking
  SAVE_FOLDER = paste0(EXPERIENCE_FOLDER, "/results/splineop/")
  if (!dir.exists(SAVE_FOLDER)) {
    dir.create(SAVE_FOLDER, recursive = TRUE)
  }
  # setup variables for file name
  Kfill = sprintf("%0*d", 2, args$K)
  samplefill = sprintf("%0*d", 2, args$sample)
  Nfill = sprintf("%0*d", 4, args$N_SAMPLES)

  # To load data
  DATA_FOLDER = paste0(dirname(EXPERIENCE_FOLDER),"/data/synth")
  #clean_signal = read.csv(paste0(DATA_FOLDER,"/raw/","K",args$K,"/K",Kfill,"id", samplefill,"N",Nfill, ".csv"),header=FALSE)
  noised_signal =  read.csv(paste0(DATA_FOLDER,"/raw/","K",args$K,"/K",Kfill,"id", samplefill,"N",Nfill,"SNR",args$snr, ".csv"),header=FALSE)

  # Set up parameters for creating SplineOP object
  noised_signal = t(noised_signal) # transpose to have the required format
  estimated_std_dev = sdHallDiff2(noised_signal)
  len_speed_estimators = as.numeric(c(20, 40, 60, 80, 100))
  n_changepoints_pred = args$K-1
  states_seed = args$sample

  # Predict pipeline
  spop <- new(SplineOP_constrained
              ,noised_signal # Data to fit
              ,len_speed_estimators # Vector fo speeds estimation
              ,estimated_std_dev # For state generation
              ,n_changepoints_pred # nb of changepoints, one less than segments
              ,states_seed # for state generation
              )
  FILE_NAME = paste0("changepoint_N",Nfill,"-K",Kfill,"-SNR",args$snr,"-ID",samplefill,".json")
  SAVE_PATH = paste0(SAVE_FOLDER/FILE_NAME)
  OUT_FILE_PATH = paste0(SAVE_PATH, FILE_NAME)
  
  results_list <- as.list(args)
  # Add the main result vector
  results_list$changepoints <- spop$get_changepoints 
  jsonlite::write_json(results_list, OUT_FILE_PATH, pretty = TRUE)
}

# 4. Execute
main(opt)
