add_noise <- function(signal,snr=30){

  sigma_signal <- sqrt(var(signal))
  noise_sd = 10^(-snr/20) * sigma_signal

  noise <- rnorm(length(signal), mean = 0, sd = noise_sd)
  noised_signal = signal + noise
  return(noised_signal)
}


