# Example of running mmsig.
# Fit previously defined reference mutational signatures to example data from a previously published cohort of multiple myeloma whole genomes
# v 09.15.2019

# Import data
muts_example <- read.delim("../data/example_muts_96.txt", stringsAsFactors = F, header=T, row.names = 1) 

muts_example_5cols <- read.delim("../data/example_muts_5cols.txt", stringsAsFactors = F, header=F) 
names(muts_example_5cols) <- c("sample", "chr", "pos", "ref", "alt")
muts_example_5cols$chr <- paste0("chr", muts_example_5cols$chr)

# Import signature reference
sig_ref <- read.delim("../data/mm_signature_definitions.txt", stringsAsFactors = F, header=T) 

# Source functions
source("util/main.R")
source("util/fitting.R")
source("util/helpers.R")
source("util/plotting.R")
source("util/strandbias.R")
source("util/bootstrap.R")

sig_out <- mm_fit_signatures(muts.input=muts_example_5cols, 
                             sig.input=sig_ref,
                             input.format = "vcf",
                             sample.sigt.profs = NULL, 
                             strandbias = TRUE,
                             bootstrap = TRUE,
                             iterations = 20,
                             refcheck=TRUE,
                             dbg=FALSE) 

plot_signatures(sig_out$estimate)

bootSigsPlot(filter(sig_out$bootstrap, sample %in% c("PD26412a", "PD26411c", "PD26414a")))

head(sig_out$strand_bias_all_3nt)

head(sig_out$strand_bias_mm1)

