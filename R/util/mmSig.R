library(seqinr)                # comp()
library(MutationalPatterns)    # cos_sim_matrix()
library(dplyr)
library(reshape2)
library(ggplot2)
library(RColorBrewer)

rotatedAxisElementText = function(angle,position='x'){
    angle     = angle[1]; 
    position  = position[1]
    positions = list(x=0,y=90,top=180,right=270)
    if(!position %in% names(positions))
        stop(sprintf("'position' must be one of [%s]",paste(names(positions),collapse=", ")),call.=FALSE)
    if(!is.numeric(angle))
        stop("'angle' must be numeric",call.=FALSE)
    rads  = (angle - positions[[ position ]])*pi/180
    hjust = 0.5*(1 - sin(rads))
    vjust = 0.5*(1 + cos(rads))
    element_text(angle=angle,vjust=vjust,hjust=hjust)
}

scale_fill_sigs <- function(...){
    ggplot2:::manual_scale('fill', 
                           values = setNames(c(RColorBrewer::brewer.pal(8, "Dark2")),
                                             c("SBS1", "SBS2", "SBS5", "SBS8", "SBS9", "SBS13", "SBS18", "SBS-MM1")), 
                           ...)
}

plot_signatures = function(sig, filename = NULL){
  "
  Plotting function for mmSig package
  filename: path to save plot. If NULL (default), returns plot to environment
  "
  if(is.null(filename)){
      par(mfrow=c(1,1), mar=c(7,5,5,10), xpd=T)
      barplot(as.matrix(t(sig[,-ncol(sig)])), 
              col = c(RColorBrewer::brewer.pal(8, "Dark2"), RColorBrewer::brewer.pal(4, "Set3")), 
              las=2, cex.names = 0.7, names.arg = sig$sampleID, space=rep(0, nrow(sig)), border = NA)
      legend("topright",legend=(colnames(sig)[-ncol(sig)]),bty="n", pch=15, 
            col=c(RColorBrewer::brewer.pal(8, "Dark2"), RColorBrewer::brewer.pal(4, "Set3")),
            cex=1, pt.cex=1, inset=c(-0.15,0.0),x.intersp = 1,y.intersp = 1)
  } else{
      pdf(filename, width = 12, height = 7)
      barplot(as.matrix(t(sig[,-ncol(sig)])), 
              col = c(RColorBrewer::brewer.pal(8, "Dark2"), RColorBrewer::brewer.pal(4, "Set3")), 
              las=2, cex.names = 0.7, names.arg = sig$sampleID, space=rep(0, nrow(sig)), border = NA)
      legend("topright",legend=(colnames(sig)[-ncol(sig)]),bty="n", pch=15, 
            col=c(RColorBrewer::brewer.pal(8, "Dark2"), RColorBrewer::brewer.pal(4, "Set3")),
            cex=1, pt.cex=1, inset=c(-0.15,0.0),x.intersp = 1,y.intersp = 1)
      dev.off()
  }
}

mm_fit_signatures = function(muts.input="../data/example_muts_96.txt", 
                             sig.input="../data/mm_signature_definitions.txt",
                             out.file="../data/mm_fitted_signatures.txt", 
                             sample.sigt.profs=NULL, 
                             dbg=FALSE) {

  "
  muts.input: samples as columns and the 96 mutational classes as rows
  if path: tab separated file
  if data frame: mutational classes as row names. Samples and total mutation count as columns. 
  sig.input: mutational signature reference with mutational classes as rows and signature exposures as columns
  if path: tab separated file
  if data frame: with the following columns: Substitution.Type, Trinucleotide, signature.1, ... signature.n
  out.file: path for output file to be generated. If NULL, returns data frame
  sample.sigt.profs: NULL = uses the hard coded signature reference for all samples. 
  Optionally provide path to list with signatures to consider for each sample, in .rds format
  dbg: FALSE = silent; TRUE = verbose
  "
  
  spit = function (dbg, mess, ...) {
    if (dbg) { print( sprintf(mess,...) ) }
  }
  
  spat = function (dbg, name, var) {
    if (dbg) {print(name); print(var)}
  }
  
  oppstrand = function(x) {   # comp() complements dna sequences
    trin1 = paste(comp(rev(strsplit(substr(x,1,3),"")[[1]]), forceToLower=F), collapse="")
    trin2 = paste(comp(rev(strsplit(substr(x,5,7),"")[[1]]), forceToLower=F), collapse="")
    return(sprintf("%s>%s",trin1,trin2))
  }
  
  # EM algorithm to estimate the signature contribution
  em_signatures = function(sigts.defn,mut.freqs,max.iter,dbg) {
    nof.sigts.defn <- ncol(sigts.defn)
    alpha = runif(nof.sigts.defn); alpha=alpha/sum(alpha) # Random start (seems to give ~identical results)
    for (i in 1:max.iter) {
      contr = t(array(alpha, dim=c(nof.sigts.defn,96))) * sigts.defn
      probs = contr/array(rowSums(contr), dim=dim(contr))
      probs[is.na(probs)] = 0
      probs = probs * mut.freqs
      old_alpha = alpha
      alpha = colSums(probs)/sum(probs)
      if (sum(abs(alpha-old_alpha))<1e-5) {
        break
      }
    }
    spit(dbg, "em: exit iteration: %d", i)
    
    return( alpha )
  }
  
  ## Read/assign input data
  
  # Mutational signature reference
  if(is.data.frame(sig.input)){
    consigts.defn <- sig.input
  } else if(file.exists(sig.input)){
    consigts.defn <- read.delim(sig.input, stringsAsFactors = F, header=T) 
  } else{
    stop("No signature reference provided, please provide either a data frame or file path")
  }
  
  # Input mutation data
  if(is.data.frame(muts.input)){
    # Input sample data: rows = 96 + 1 3nuc mutation patterns (+ total), cols = 89 samples + 1 total
    samples.muts <- muts.input
  } else if(file.exists(muts.input)){
    samples.muts <- read.delim(muts.input, stringsAsFactors = F, header=T, row.names=1) 
  } else{
    stop("No input mutation data provided, please provide either a data frame or file path")
  }
  
  # Process signature reference
  tcons.defn <- t(consigts.defn[,3:ncol(consigts.defn)])
  mutlist = paste(consigts.defn[,2],paste(substr(consigts.defn[,2],1,1),substr(consigts.defn[,1],3,3),substr(consigts.defn[,2],3,3),sep=""),sep=">")
  # muttype = rep(1:96,2)
  # mutlist_oppstrand = sapply(mutlist, oppstrand)
  # muttype = rep(1:96,2)
  # names(muttype) = c(mutlist,mutlist_oppstrand)
  consigts.defn <- sapply(consigts.defn[,3:ncol(consigts.defn)], as.numeric)  #n be cautious in the original the first two columns are included
  consigt.names <- colnames(consigts.defn)
  rownames(consigts.defn)<-mutlist
  
  # Process sample mutational profiles
  samples.muts <- samples.muts[names(samples.muts) != "Total"]                       # remove the totals column
  samples.muts <- samples.muts[1:96,]                                                # remove additional rows (e.g. "Total)"
  samples <- colnames(samples.muts)

  # Assign which signatures to fit to each sample
  if (isTRUE(file.exists(file.path(sample.sigt.profs)))) {                           # workaround to prevent crash if NULL
      library("R.filesets")            # loadRDS()
      spit(dbg, "using mm signature profiles from file: %s", sample.sigt.profs)
      sigt.profs <- loadRDS(sample.sigt.profs)
      
  } else if (is.list(sample.sigt.profs)){
      spit(dbg, "using mm signature profiles from input argument")
      sigt.profs <- sample.sigt.profs
          
  } else {
      spit(dbg, "using prior multiple myeloma signature profiles")
      mm.sigts <- (c("SBS1","SBS2","SBS5","SBS8","SBS9",
                     "SBS13","SBS18","SBS.MM1"))
      
      sigt.profs <- list()
      for (i in 1:length(samples)) {
          sigt.profs[[samples[i]]] <- mm.sigts
      }
  }
  
  # Mutational signature fitting procedure for each infividual sample
  sigt.fraction = array(NA,dim=c(length(consigt.names), length(samples)))
  rownames(sigt.fraction) <- consigt.names
  colnames(sigt.fraction) <- samples
  max.em.iter=2000
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
    sample.cos.sim.meas <- cos_sim_matrix(reconstructed, matrix(sample.mut.freqs, ncol=1))
    spat(dbg, "sample.cos.sim.meas", sample.cos.sim.meas)
    
    rem.alpha <- sampleAlpha                     # holds the final result
    rem.sample.consigts.defn <- sample.consigts.defn
    spit(dbg, "length of rem.alpha: %d", length(rem.alpha))
    reducing = TRUE
    
    # Signature profile shrinkage by cosine similarity (removing signatures that are not necessary to explain profile)
    while (reducing) { 
      spat(dbg, "in the while, rem.alpha: ", rem.alpha)
      cosReduction <- NULL
      rem.names <- setdiff(names(rem.alpha),c("SBS1","SBS5"))
      if(length(rem.names) == 0){ ## Avoiding script crash when only SBS1 and 5 are present.
        spit(dbg, "removed all signatures except SBS1 and SBS5: exiting while...")
        break
      }
      for(c in rem.names){
        spit(dbg, "doing c: %s", c)
        red.sample.consigts.defn <- rem.sample.consigts.defn[,colnames(rem.sample.consigts.defn)!=c]
        red.alpha <- em_signatures(sigts.defn=red.sample.consigts.defn,mut.freqs=sample.mut.freqs,max.iter=max.em.iter,dbg=dbg)
        red.reconstructed <- red.sample.consigts.defn %*% red.alpha * sum(sample.mut.freqs)
        red.cos.sim.meas <- cos_sim_matrix(red.reconstructed, matrix(sample.mut.freqs, ncol=1))
        cosReduction <- c(cosReduction, sample.cos.sim.meas-red.cos.sim.meas)
      }
      names(cosReduction) <- rem.names
      if (min(cosReduction) < 0.01) {
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

  if(is.null(out.file)){
    return(sig)
  } else{
    write.table(sig, out.file, sep="\t", quote=F)
    figname <- sub("\\.[a-z]{3}$", "_plot.pdf", out.file)
    plot_signatures(sig, filename = figname)
  }
}



bootstrap_mm_signatures <- function(muts.input, 
                                    sig.input,
                                    sample.sigt.profs=NULL, 
                                    iterations = 1000){
    "
    Bootstrapping function for mm_fit_signatures
    Draw b = iterations mutational profiles for each sample from its multinomial distribution
    Performs signature fitting independently for each mutational profile
    For each mutational signature, returns its relative contribution as point estimate and bootstrapping mean with 95 % CI
    Can take a sample.sigt.profs argument that is passed directly to mm_fit_signatures. 
    "
    # Generate point estimates
    sig_est <- mm_fit_signatures(muts.input=muts.input, 
                                 sig.input=sig.input,
                                 out.file=NULL,
                                 sample.sigt.profs=sample.sigt.profs, 
                                 dbg=FALSE)
    
    sig_est$sample <- row.names(sig_est)
    sig_est <- melt(sig_est[names(sig_est) != "mutations"], 
                    id.vars = "sample", 
                    variable.name = 'signature', 
                    value.name = 'estimate')
    
    sig_est$signature <- as.character(sig_est$signature)
    
    # Setup
    samples <- names(muts.input)[names(muts.input) != 'Total']
    classes <- row.names(muts.input)[row.names(muts.input) != 'Total']
    
    # List to populate with signatures
    mutSigs <- list()
    
    for(i in 1:length(samples)){
        # Loop through samples, generating a data frame of signature contributions for each
        sub <- as.integer(muts.input[classes,i])
        total <- sum(sub)
        
        # sample new 96-classes profiles from the multinomial distribution
        bootMat <- data.frame(rmultinom(iterations, total, sub/total))
        row.names(bootMat) <- classes
        
        # prepare the signatures to fit for each sample
        if(is.list(sample.sigt.profs)){
            sigt.prof <- list()
            for(s in 1:ncol(bootMat)){
                sigt.prof[[names(bootMat)[s]]] <- sample.sigt.profs[[samples[i]]]
            }
        } else {
            sigt.prof <- NULL
        }
        
        ### Run mmSig
        sig_out <- mm_fit_signatures(muts.input=bootMat, 
                                     sig.input=sig.input,
                                     out.file = NULL,
                                     sample.sigt.profs=sigt.prof, 
                                     dbg=FALSE)
        
        mutSigs[[i]] <- sig_out
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
    mutSigsSummary <- left_join(mutSigsSummary, sig_est, by = c("sample", "signature"))
    
    return(mutSigsSummary)
}


bootSigsPlot <- function(mutSigsSummary){
    ggplot(mutSigsSummary, aes(signature, estimate, fill = signature))+
        geom_bar(position="dodge", stat="identity")+
        geom_errorbar(data = mutSigsSummary, mapping = aes(x = signature, ymin = CI025, ymax = CI975))+
        facet_grid(~sample)+
        theme(text = element_text(size = 12),
              axis.text.x = rotatedAxisElementText(90, 'top'),
              axis.title.x = element_blank(),
              strip.background = element_blank(),
              legend.position = 'none')+
        scale_fill_sigs()+
        labs(y = 'Relative contribution')
}
