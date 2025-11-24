library(splineOP)


SAVE_DIRECTORY = dirname(rstudioapi::getSourceEditorContext()$path)
LOAD_DIRECTORY = paste0(SAVE_DIRECTORY, "/K", K_segments,"/raw/")
SNR = 30

for (K_segments in c(5,10, 15))
  {#
  LOAD_DIRECTORY = paste0(SAVE_DIRECTORY, "/K", K_segments,"/raw/")
  NOISE_DIRECTORY = paste0(SAVE_DIRECTORY, "/K", K_segments,"/noised/")
  if (!dir.exists(NOISE_DIRECTORY)) {
    dir.create(NOISE_DIRECTORY, recursive = TRUE, showWarnings = FALSE)
  }
  for (sample in 1:100)
    {
    # define file path
    Kfill = sprintf("%0*d", 2, K_segments)
    samplefill = sprintf("%0*d", 2, sample)
    Nfill = sprintf("%0*d", 4, N)
    file_path = paste0(LOAD_DIRECTORY,"K",Kfill,"id", samplefill,"N",Nfill, ".csv")
    signal = read.csv(file = file_path, header=FALSE)    # Load data

    # add noise
    noised_signal <- add_noise(signal, snr=SNR)
    # save noised file
    write.table(noised_signal,
                paste0(NOISE_DIRECTORY,"K",Kfill,"id", samplefill,"N",Nfill,"SNR",SNR,".csv"),
                sep = ",",
                col.names=FALSE,
                row.names=FALSE)
    }
  }
