#' Clean input VCF-like file
#'
#' @param snvs input mutation catalog in 5-columns format
#'
#' @return clean mutation catalog in 5-columns format
#'
vcfCleaner <- function(snvs){

  snvs <- snvs[snvs$ref %in% c("A","C","G","T") & snvs$alt %in% c("A","C","G","T"),] # remove indels
  snvs$chr <- as.character(snvs$chr)
  snvs <- snvs[snvs$chr %in% paste0("chr", c(1:22, "X")),]
  snvs$pos <- as.numeric(snvs$pos)
  snvs <- snvs[stats::complete.cases(snvs),]
  return(snvs)
}
