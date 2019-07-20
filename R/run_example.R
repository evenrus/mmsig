# Running mmSig: fit a previously defined reference mutational signatures to example data from a previously published cohort of multiple myeloma whole genomes
# The mmSig package has options to run with input/output from files. Here, we run the package as it would be included in a larger analysis script. 
# v 07.20.2019

# Import data
muts_example <- read.delim("../data/example_muts_96.txt", stringsAsFactors = F, header=T, row.names=1) 

# Import signature reference
sig_ref <- read.delim("../data/mm_signature_definitions.txt", stringsAsFactors = F, header=T, row.names=1) 

# Source functions
source("util/mmSig.R")

sig_out <- mm_fit_signatures(muts.input=muts_example, 
                             sig.input=sig_ref,
                             out.file=NULL, 
                             sample.sigt.profs=NULL, 
                             dbg=FALSE,
                             return.fig = TRUE) 

plot_signatures(sig_out)