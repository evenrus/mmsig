#' Bootstrapping function for mmsig
#'
#' @param samples.muts 96 classes mutational profile for all samples
#' @param consigts.defn 96 classes mutational profile for mutational signature reference
#' @param sigt.profs which signatures to fit for each sample
#' @param cos_sim_threshold cosine similarity threshold to remove signatures
#' @param force_include vector with the names of signatures to always keep in the final profile of every sample
#' @param iterations number of mutational profiles to draw from the multinomial distribution of each sample
#' @importFrom plyr create_progress_bar
#' @importFrom plyr progress_text
#' @importFrom dplyr bind_rows
#'
#' @return For each mutational signature, returns its relative contribution as point estimate and bootstrapping mean with 95 % CI
#'
bootstrap_fit_signatures <- function(samples.muts,
                                     consigts.defn,
                                     sigt.profs,
                                     cos_sim_threshold,
                                     force_include,
                                     iterations = 1000){

  # Setup
  samples <- names(samples.muts)
  classes <- row.names(samples.muts)

  # List to populate with signatures
  mutSigs <- list()

  # Setup progress bar
  pbar <- create_progress_bar('text')
  pbar$init(length(samples))

  for(i in 1:length(samples)){
    # Loop through samples, generating a data frame of signature contributions for each
    sub <- as.integer(samples.muts[classes,i])
    total <- sum(sub)
    # sample new 96-classes profiles from the multinomial distribution
    bootMat <- data.frame(rmultinom(n = iterations, size = total, prob = sub/total))
    row.names(bootMat) <- classes

    # prepare the signatures to fit for each sample
    sig.prof <- list()

    for(s in 1:ncol(bootMat)){
      sig.prof[[names(bootMat)[s]]] <- sigt.profs[[samples[i]]]
    }

    ### Run mmSig
    sig_out <- fit_signatures(samples.muts=bootMat,
                              consigts.defn=consigts.defn,
                              sigt.profs=sig.prof,
                              cos_sim_threshold=cos_sim_threshold,
                              force_include=force_include,
                              dbg=FALSE)

    mutSigs[[i]] <- sig_out
    pbar$step()
  }

  names(mutSigs) <- samples

  # Generate final summary data frame

  ## Summary statistics
  my_summary <- function(x){
    c(mean(x), quantile(x, probs = 0.025), quantile(x, probs = 0.975))
  }

  mutSigsSummary <- list()
  for(i in 1:length(mutSigs)){
    s <- names(mutSigs)[i]
    temp <- mutSigs[[i]]
    out <- data.frame(t(sapply(temp[names(temp) != "mutations"], my_summary)))
    names(out) <- c('mean', 'CI025', 'CI975')
    out$signature <- row.names(out)
    out$sample <- s
    out <- out[c('sample', 'signature', 'mean', 'CI025', 'CI975')]
    mutSigsSummary[[i]] <- out
  }

  mutSigsSummary <- bind_rows(mutSigsSummary)

  return(mutSigsSummary)
}
