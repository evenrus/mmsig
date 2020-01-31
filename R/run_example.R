# Example of running mmsig.
# Fit previously defined reference mutational signatures to example data from a previously published cohort of multiple myeloma whole genomes
# v 01.03.2020

# Import data
# mmsig can take input data both as a matrix of 96 mutational classes profiles, and in a 5-column format where each row is a mutation with chromosome, position, etc. 
# strand bias can only be analyzed from 5-column input
# the current version of mmsig only supports 5-column format data aligned to hg19/GRCh37

muts_example <- read.delim("../data/example_muts_96.txt", stringsAsFactors = F, header=T, row.names = 1) 
muts_example_5cols <- read.delim("../data/example_muts_5cols.txt", stringsAsFactors = F, header=T) 

# Import signature reference
sig_ref <- read.delim("../data/mm_signature_definitions.txt", stringsAsFactors = F, header=T) 

# Source functions
source("util/main.R")
source("util/fitting.R")
source("util/helpers.R")
source("util/plotting.R")
source("util/strandbias.R")
source("util/bootstrap.R")

# fit signatures -- runs in < 5 minutes on a MacBook pro with 2.8 GHz Intel Core i7 processor.
# Bootstrapping large datasets with many iterations can significantly increase runtime. 

sig_out <- mm_fit_signatures(muts.input=muts_example_5cols, 
                             sig.input=sig_ref,
                             input.format = "vcf",
                             sample.sigt.profs = NULL, 
                             strandbias = TRUE,
                             bootstrap = TRUE,
                             iterations = 20, # 1000 iterations recommended for stable results
                             refcheck=TRUE,
                             dbg=FALSE) 

# data summaries and plots
plot_signatures(sig_out$estimate)

bootSigsPlot(filter(sig_out$bootstrap, sample %in% c("PD26412a", "PD26411c", "PD26414a")))

head(sig_out$strand_bias_all_3nt)

head(sig_out$strand_bias_mm1)

head(sig_out$strand_bias_SBS35)



