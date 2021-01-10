#' Mutational signature fitting with mmsig
#'
#' @param muts.input mutation input data frame as specified under input.format
#' @param sig.input mutational signature reference with mutational classes as rows and signature exposures as columns: Substitution.Type, Trinucleotide, signature.1, ... signature.n
#' @param input.format vcf: five column vcf-like data frame with the following columns: sample, chr (e.g. chr1), pos, ref, alt. classes: samples as columns and the 96 mutational classes as rows
#' @param sample.sigt.profs NULL = use all signatures provided in the reference. Optionally provide list with signatures to consider for each sample
#' @param bootstrap TRUE/FALSE for whether bootstrapping is to be performed
#' @param iterations number of bootstrapping iterations to perform (only if bootstrap == TRUE)
#' @param strandbias TRUE/FALSE for whether transcriptional strand bias should be tested for (only for vcf-like input format)
#' @param refcheck check that input mutational catalog (if vcf-format) is aligned to hg19
#' @param cos_sim_threshold cosine similarity threshold below which signatures are removed from the final profile
#' @param force_include vector with the names of signatures to always keep in the final profile of every sample
#' @param dbg FALSE = silent; TRUE = verbose
#' @importFrom dplyr left_join
#' @importFrom deconstructSigs mut.to.sigs.input
#' @import BSgenome.Hsapiens.UCSC.hg19
#'
#' @return mutational signature fitting results for all samples
#' @export
#'
mm_fit_signatures = function(muts.input,
                             sig.input,
                             input.format = "vcf",
                             sample.sigt.profs=NULL,
                             bootstrap=FALSE,
                             iterations=1000,
                             strandbias=FALSE,
                             refcheck=TRUE,
                             cos_sim_threshold=0.01,
                             force_include=c("SBS1", "SBS5"),
                             dbg=FALSE) {
  "

    "

  options(scipen = 999)

  #####################################################
  ########        Read/assign input data       ########
  #####################################################

  # Mutational signature reference
  consigts.defn <- sig.input

  # Input mutation data
  if(input.format == "vcf") {

    if(sum(c("sample", "chr", "pos", "ref", "alt") %in% names(muts.input)) != 5){
      stop("ERROR: Please provide the following data columns: sample, chr, pos, ref, alt")
    }

    muts.input <- vcfCleaner(muts.input)

    if(nrow(muts.input) == 0){
      stop("ERROR: Inappropriate input data format")
    }

    if(refcheck){
      if(!refCheck(muts.input, BSgenome.Hsapiens.UCSC.hg19)){
        stop("ERROR: Wrong reference genome, please provide mutational data aligned to hg19/GRCh37")
      }
    }

    # Input sample data
    samples.muts <- mut.to.sigs.input(mut.ref = muts.input,
                                      sample.id = "sample",
                                      chr = "chr",
                                      pos = "pos",
                                      ref = "ref",
                                      alt = "alt",
                                      bsg = BSgenome.Hsapiens.UCSC.hg19)

    samples.muts <- as.data.frame(t(samples.muts))

  } else if(input.format == "classes"){
    samples.muts <- muts.input
    samples.muts <- samples.muts[names(samples.muts) != "Total"]                       # remove totals column
    samples.muts <- samples.muts[1:96,,drop=FALSE]                                     # remove additional rows (e.g. "Total)"
  } else{
    stop("ERROR: Invalid input format, please provide a vcf-like format ('vcf') or 96 mutational classes ('classes')")
  }

  # Process signature reference
  tcons.defn <- t(consigts.defn[,3:ncol(consigts.defn)])
  mutlist = paste(consigts.defn[,2],paste(substr(consigts.defn[,2],1,1),substr(consigts.defn[,1],3,3),substr(consigts.defn[,2],3,3),sep=""),sep=">")
  consigts.defn <- sapply(consigts.defn[,3:ncol(consigts.defn)], as.numeric)  #n be cautious in the original the first two columns are included
  rownames(consigts.defn)<-mutlist
  ref_signatures <- colnames(consigts.defn) # names of signatures in reference

  # Process sample mutational profiles
  samples <- colnames(samples.muts)

  # Assign which signatures to fit to each sample
  if(is.list(sample.sigt.profs)){
    spit(dbg, "using mm signature profiles from input argument")
    sigt.profs <- sample.sigt.profs
  } else {
    spit(dbg, "defaulting to use all signatures in the provided reference")
    mm.sigts <- ref_signatures

    sigt.profs <- list()
    for (i in 1:length(samples)) {
      sigt.profs[[samples[i]]] <- mm.sigts
    }
  }

  # Define output variable
  output <-list()

  #####################################################
  ########       Fit mutational signatures     ########
  #####################################################

  sigfit <- fit_signatures(samples.muts=samples.muts,
                           consigts.defn=consigts.defn,
                           sigt.profs=sigt.profs,
                           cos_sim_threshold=cos_sim_threshold,
                           force_include=force_include,
                           dbg=dbg)

  output$estimate <- sigfit

  #####################################################
  ########             Bootstrapping           ########
  #####################################################

  if(bootstrap){

    # Point estimates
    sig_est <- sigfit
    sig_est$sample <- row.names(sig_est)
    sig_est <- melt(sig_est[names(sig_est) != "mutations"],
                    id.vars = "sample",
                    variable.name = 'signature',
                    value.name = 'estimate')

    sig_est$signature <- as.character(sig_est$signature)

    sigboot <- bootstrap_fit_signatures(samples.muts = samples.muts,
                                        consigts.defn = consigts.defn,
                                        sigt.profs = sigt.profs,
                                        iterations = iterations,
                                        cos_sim_threshold = cos_sim_threshold,
                                        force_include = force_include)

    sigboot <- left_join(sigboot, sig_est, by = c("sample", "signature"))

    output$bootstrap <- sigboot
  }

  #####################################################
  ########              Strand Bias            ########
  #####################################################

  if(strandbias & input.format == "vcf"){

    strand_bias_out <- getStrandBias(muts.input)

    output$strand_bias_all_3nt <- strand_bias_out$all_3nt
    output$strand_bias_mm1 <- strand_bias_out$mm1
    output$strand_bias_SBS35 <- strand_bias_out$SBS35

  } else if(strandbias & input.format != "vcf") {
    warning("Transcriptional strand bias cannot be estimated from 96 classes input. \n Please provide vcf-like input format.")
  }

  output$mutmatrix <- samples.muts

  return(output)
}
