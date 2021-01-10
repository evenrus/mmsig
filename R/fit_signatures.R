#' Signature fitting function for mmsig
#'
#' @param samples.muts 96 classes mutational profile for all samples
#' @param consigts.defn 96 classes mutational profile for mutational signature reference
#' @param sigt.profs which signatures to fit for each sample
#' @param cos_sim_threshold cosine similarity threshold to remove signatures
#' @param force_include vector with the names of signatures to always keep in the final profile of every sample
#' @param dbg boolean whether to print or not
#'
#' @return estimated mutational signature profile of all samples
#'
fit_signatures = function(samples.muts,
                          consigts.defn,
                          sigt.profs,
                          cos_sim_threshold,
                          force_include,
                          dbg=dbg) {
  max.em.iter=2000
  consigt.names <- colnames(consigts.defn)
  samples <- colnames(samples.muts)

  # Mutational signature fitting procedure for each individual sample
  sigt.fraction = array(NA,dim=c(length(consigt.names), length(samples)))
  rownames(sigt.fraction) <- consigt.names
  colnames(sigt.fraction) <- samples

  for (j in 1:length(samples)) {
    sample.mut.freqs = as.numeric(samples.muts[,j])
    sample.mut.freqs[is.na(sample.mut.freqs)] = 0
    sample.sigts <- unique(sigt.profs[[ samples[j] ]])
    sample.sigts <- sample.sigts[match(consigt.names[consigt.names %in% sample.sigts], sample.sigts)]
    sample.consigts.defn  <- consigts.defn[, colnames(consigts.defn) %in% sample.sigts]
    spat(dbg, "colnames sample.consigts.defn (before em)", colnames(sample.consigts.defn))
    alpha <- em_signatures(sigts.defn=sample.consigts.defn,mut.freqs=sample.mut.freqs,max.iter=max.em.iter,dbg=dbg)
    spat(dbg, "alpha", alpha)
    sample.consigts.defn <- sample.consigts.defn[, colnames(sample.consigts.defn) %in% names(alpha)]   # output sample.... should be identical to input sample....
    sampleAlpha <- alpha[match(colnames(sample.consigts.defn), names(alpha))]

    if (!all(alpha==sampleAlpha)) {stop("non-identical alphas")}

    spat(dbg, "colnames: sample.consigts.defn (after em (and reduction))", colnames(sample.consigts.defn))
    reconstructed <- sample.consigts.defn %*% alpha * sum(sample.mut.freqs)
    sample.cos.sim.meas <- MutationalPatterns::cos_sim_matrix(reconstructed, matrix(sample.mut.freqs, ncol=1))
    spat(dbg, "sample.cos.sim.meas", sample.cos.sim.meas)

    rem.alpha <- sampleAlpha                     # holds the final result
    rem.sample.consigts.defn <- sample.consigts.defn
    spit(dbg, "length of rem.alpha: %d", length(rem.alpha))
    reducing = TRUE

    # Signature profile shrinkage by cosine similarity (removing signatures that are not necessary to explain profile)
    while (reducing) {
      spat(dbg, "in the while, rem.alpha: ", rem.alpha)
      cosReduction <- NULL
      rem.names <- setdiff(names(rem.alpha), force_include)
      if(length(rem.names) == 0){ ## Avoiding script crash when only the forced signatures are present.
        spit(dbg, "removed all signatures except forced inclusion: exiting while...")
        break
      }

      for(c in rem.names){
        spit(dbg, "doing c: %s", c)
        red.sample.consigts.defn <- rem.sample.consigts.defn[,colnames(rem.sample.consigts.defn)!=c]
        red.alpha <- em_signatures(sigts.defn=red.sample.consigts.defn,mut.freqs=sample.mut.freqs,max.iter=max.em.iter,dbg=dbg)
        red.reconstructed <- red.sample.consigts.defn %*% red.alpha * sum(sample.mut.freqs)
        red.cos.sim.meas <- MutationalPatterns::cos_sim_matrix(red.reconstructed, matrix(sample.mut.freqs, ncol=1))
        cosReduction <- c(cosReduction, sample.cos.sim.meas-red.cos.sim.meas)
      }
      names(cosReduction) <- rem.names
      if (min(cosReduction) < cos_sim_threshold) {
        spit(dbg, "removing: %s", names(cosReduction)[which.min(cosReduction)])
        rem.sample.consigts.defn <- rem.sample.consigts.defn[,- which(colnames(rem.sample.consigts.defn)==names(which.min(cosReduction)))]
        rem.alpha <-  em_signatures(sigts.defn=rem.sample.consigts.defn,mut.freqs=sample.mut.freqs,max.iter=max.em.iter,dbg=dbg)
        reducing = TRUE
      }
      else {
        spit(dbg, "exiting while...")
        reducing = FALSE
      }
    }

    spit(dbg,"... while exited")
    rem.alpha.names <- names(rem.alpha)

    for (n in 1:length(consigt.names)) {
      if (consigt.names[n] %in% rem.alpha.names) {
        sigt.fraction[n,j] <- rem.alpha[consigt.names[n]]
      }
      else {
        sigt.fraction[n,j] <- 0
      }
    }
  }
  spat(dbg, "sigt.fraction", sigt.fraction)
  tdf.sigt.fraction <- as.data.frame(t(sigt.fraction))
  colsums.samples.muts <- colSums(samples.muts)
  sig <- cbind(tdf.sigt.fraction, "mutations"=colsums.samples.muts)
  colnames(sig)[colnames(sig)=="SBS.MM1"]<-"SBS-MM1"

  return(sig)
}

