# mmsig
A method for fitting a known mutational signature reference to mutational catalogues from cancer samples

Accompanying the publication: 
Timing the initiation of multiple myeloma https://papers.ssrn.com/sol3/papers.cfm?abstract_id=3409453

Authors: Even Holth Rustad, Francesco Maura, Nicos Angelopoulos and Venkata Yellapantula

Usage examples can be found in ./R/run_example.R, applying mutational signature analysis to a set of multiple myeloma whole genome sequences. 

Developed and tested in R version 3.6.1. using the following package versions:
seqinr_3.6-1
MutationalPatterns_1.10.0
plyr_1.8.4
dplyr_0.8.3
reshape2_1.4.3
ggplot2_3.2.1
RColorBrewer_1.1-2
tidyr_1.0.0
deconstructSigs_1.8.0
BSgenome_1.52.0
TxDb.Hsapiens.UCSC.hg19.knownGene_3.2.2
BSgenome.Hsapiens.UCSC.hg19_1.4.0 