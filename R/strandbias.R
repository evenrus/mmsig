#' Statistical tests for transcriptional strand bias
#'
#' @param data_5cols input mutational catalog in 5-column format
#' @return list of results from strand bias in each trinucleotide context and pooled analysis for SBS-MM1 and SBS35
#' @importFrom dplyr %>%
#' @importFrom tibble rownames_to_column
#' @importFrom reshape2 melt
#' @importFrom tidyr separate
#' @importFrom dplyr filter
#' @importFrom MutationalPatterns strand_bias_test
#' @importFrom reshape2 dcast
#' @importFrom dplyr rowwise
#' @importFrom stats poisson.test
#' @importFrom TxDb.Hsapiens.UCSC.hg19.knownGene TxDb.Hsapiens.UCSC.hg19.knownGene
#' @importFrom GenomicRanges GRanges
#' @importFrom IRanges IRanges
#' @import BSgenome.Hsapiens.UCSC.hg19
#'

getStrandBias <- function(data_5cols){

  # Generate GRange list from data in vcf-like format
  sample_list <- list()
  samples <- unique(data_5cols$sample)

  for(i in (1:length(samples))) {
    temp <- data_5cols[data_5cols$sample== samples[i],]
    temp_gr <- with(temp,
                    GRanges(chr,
                            IRanges(start=pos, end=pos),
                            REF= ref,
                            ALT=alt))
    genome(temp_gr) <- "hg19"

    sample_list[[i]]<- (temp_gr)
  }

  names(sample_list) <- samples

  # transcriptional strand annotation
  genes_hg19 <- suppressMessages(GenomicFeatures::genes(TxDb.Hsapiens.UCSC.hg19.knownGene))

  mut_mat_stranded <- MutationalPatterns::mut_matrix_stranded(vcf_list = sample_list,
                                                              ref_genome = "BSgenome.Hsapiens.UCSC.hg19",
                                                              ranges = genes_hg19,
                                                              mode = "transcription")

  mut_df_stranded <- as.data.frame(mut_mat_stranded) %>%
    rownames_to_column(var = "type") %>%
    melt(id.var = "type", variable.name = 'group', value.name = 'no_mutations') %>%
    separate(col = 'type', into = c('type', 'strand'), sep = "-")

  # Test for all transcriptional strand bias by 3-nt and sample
  mut_strandtest <- suppressMessages(strand_bias_test(mut_df_stranded))

  # Test for mm1
  mm1 <- c("C[C>T]A","G[C>T]A", "G[C>T]C", "G[C>T]G", "G[C>T]T")

  mut_mm1 <- mut_df_stranded %>%
    filter(type %in% mm1) %>%
    dplyr::group_by(group, strand) %>%
    dplyr::summarise(count = sum(no_mutations)) %>%
    as.data.frame() %>%
    dcast(group ~ strand, value.var = "count") %>%
    rowwise() %>%
    dplyr::mutate(ratio = transcribed/untranscribed) %>%
    as.data.frame()

  # Applying the same poisson test as the strand_bias_test function

  p_poisson <- list()
  for(i in 1:nrow(mut_mm1)){
    transcribed <- mut_mm1[i,"transcribed"]
    untranscribed <- mut_mm1[i,"untranscribed"]
    pval <- poisson.test(c(transcribed, untranscribed), r = 1)$p.value

    p_poisson[[i]] <- ifelse(is.numeric(pval), pval, NA)
  }

  mut_mm1$p_poisson <- unlist(p_poisson)

  mut_mm1 <- dplyr::mutate(mut_mm1, MM1_flag = ifelse(p_poisson < 0.05 & ratio > 1, "*", ""))

  # Test for platinum / SBS35

  SBS35 <- c("C[C>A]C", "C[C>T]C")

  mut_SBS35 <- mut_df_stranded %>%
    filter(type %in% SBS35) %>%
    dplyr::group_by(group, strand) %>%
    dplyr::summarise(count = sum(no_mutations)) %>%
    as.data.frame() %>%
    dcast(group ~ strand, value.var = "count") %>%
    rowwise() %>%
    dplyr::mutate(ratio = transcribed/untranscribed) %>%
    as.data.frame()

  # Applying the same poisson test as the strand_bias_test function

  p_poisson <- list()
  for(i in 1:nrow(mut_SBS35)){
    transcribed <- mut_SBS35[i,"transcribed"]
    untranscribed <- mut_SBS35[i,"untranscribed"]
    pval <- poisson.test(c(transcribed, untranscribed), r = 1)$p.value

    p_poisson[[i]] <- ifelse(is.numeric(pval), pval, NA)
  }

  mut_SBS35$p_poisson <- unlist(p_poisson)

  mut_SBS35 <- dplyr::mutate(mut_SBS35, SBS35_flag = ifelse(p_poisson < 0.05 & ratio > 1, "*", ""))

  # output

  output <- list(all_3nt = mut_strandtest,
                 mm1 = mut_mm1,
                 SBS35 = mut_SBS35)

  return(output)
}

