## Bootstrapping function for mmsig

bootstrap_fit_signatures <- function(samples.muts, 
                                     consigts.defn,
                                     sigt.profs=sigt.profs, 
                                     iterations = 1000){
    
    "
    Bootstrapping function for mm_fit_signatures
    Draw b = iterations mutational profiles for each sample from its multinomial distribution
    Performs signature fitting independently for each mutational profile
    For each mutational signature, returns its relative contribution as point estimate and bootstrapping mean with 95 % CI
    Can take a sample.sigt.profs argument that is passed directly to mm_fit_signatures. 
    "
    
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