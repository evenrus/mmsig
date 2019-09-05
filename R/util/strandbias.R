# Strand bias function for mmsig

## STRAND BIAS

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
        sample_list[[i]]<- (temp_gr)
    }
    
    names(sample_list) <- samples
    
    # transcriptional strand annotation
    genes_hg19 <- genes(TxDb.Hsapiens.UCSC.hg19.knownGene)
    
    mut_mat_stranded <- mut_matrix_stranded(vcf_list = sample_list, ref_genome = ref_genome, ranges = genes_hg19, mode = "transcription")
    
    mut_df_stranded <- as.data.frame(mut_mat_stranded) %>%
        dplyr::mutate(type = rownames(as.data.frame(mut_mat_stranded))) %>%
        melt(variable.name = 'group', value.name = 'no_mutations') %>%
        separate(col = 'type', into = c('type', 'strand'), sep = "-")
    
    # Test for all transcriptional strand bias by 3-nt and sample
    mut_strandtest <- strand_bias_test(mut_df_stranded)
    mut_strandtest$FDR <- p.adjust(mut_strandtest$p_poisson, method = "fdr")
    
    mut_strandtest <- dplyr::mutate(mut_strandtest, significant = ifelse(FDR < 0.1, "*", ""))
    
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
        pval <- poisson.test(transcribed, untranscribed, r = 1)$p.value
        
        p_poisson[[i]] <- ifelse(is.numeric(pval), pval, NA)
    }
    
    mut_mm1$p_poisson <- unlist(p_poisson)
    mut_mm1$FDR <- p.adjust(mut_mm1$p_poisson, method = "fdr")
    mut_mm1 <- dplyr::mutate(mut_mm1, MM1_flag = ifelse(FDR < 0.1 & ratio > 1, "*", ""))
    
    
    output <- list(all_3nt = mut_strandtest,
                   mm1 = mut_mm1)
    return(output)
}




