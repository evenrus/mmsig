# Running mmSig with bootstrapping to generate 95 % confidence intervals for mutational signature contributions
# v 08.06.2019

# Import data
muts_example <- read.delim("../data/example_muts_96.txt", stringsAsFactors = F, header=T, row.names = 1) 

# Import signature reference
sig_ref <- read.delim("../data/mm_signature_definitions.txt", stringsAsFactors = F, header=T) 

# Source functions
source("util/mmSig.R")

# Subset test dataset
muts.input <- muts_example[,1:5]

# Generates a list of which mutational signatures to fit for each sample
sig.prof <- list()
for(s in 1:ncol(muts.input)){
    sigt.prof[[names(muts.input)[s]]] <- c("SBS1", "SBS5", "SBS13")
}

# Run bootstrapping and signature fitting
mutSigsSummary <- bootstrap_mm_signatures(muts.input = muts.input, 
                                          sig.input = sig_ref,
                                          sample.sigt.profs=sig.prof, 
                                          iterations = 100)

## Plotting
bootSigsPlot(mutSigsSummary)

ggsave("../data/bootstrapped_relative_sigs_contributions.pdf", height = 5, width = 7)
