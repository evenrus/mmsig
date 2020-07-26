#' Check if the input mutation catalog is aligned to the correct reference genome
#'
#' @param snvs input mutation catalog in 5-columns format
#' @param genome BSgenome object with the full reference genome used for analysis
#'
#' @return reference genome match: TRUE or FALSE
#'
refCheck <- function(snvs, genome){
  # number of mutations to check
  nCheck <- ifelse(nrow(snvs) > 100, 100, nrow(snvs))
  # subset mutations to check
  sub <- snvs[sample(nrow(snvs), nCheck),]
  ## Get the reference genome positions
  newRef <- BSgenome::getSeq(genome, as.character(sub$chr), sub$pos, sub$pos)
  # Check reference genome positions
  out <- all(newRef == sub$ref)
  return(out)
}
