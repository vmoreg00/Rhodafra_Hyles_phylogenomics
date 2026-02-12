#!/usr/bin/env Rscript

# test if there is at least four argument: if not, return an error
args = commandArgs(trailingOnly=TRUE)
if (length(args) != 1) {
  cat("Input directory must be supplied!\n")
  cat("\tRscript vanilla 07.01_plot_twisst.R <input_dir>")
  stop("", call.=FALSE)
}
datapath <- args[1]

# check file and directory
if(!dir.exists(datapath)){
  stop("Input directory not found!")
}

# Load plot_twisst.R functions
source("~/hDNAhyles/src/twisst/plot_twisst.R")



# Get the file names
weights_file <- paste0(datapath, "/02_raxml_sw_phylo.weights.tsv.gz")
window_data_file <- paste0(datapath, "/01_raxml_sw_phylo.data.tsv")
# Import data
twisst_data <- import.twisst(weights_files = weights_file,
                             window_data_files = window_data_file)
# Smooth
twisst_data_smooth <- smooth.twisst(twisst_data, span_bp = 500000)
# Plot weights in windows
pdf(paste0(datapath, "/weights_plot.pdf"), width = 15, height = 15)
plot.twisst(twisst_data_smooth, mode = 3,
            ncol_topos = min(5, length(twisst_data$topos)), ncol_weights = 4)
dev.off()
# Plot summary barplot
pdf(paste0(datapath, "/summary.pdf"), width = 5+length(twisst_data$topos), height = 4)
plot.twisst.summary(twisst_data)
dev.off()
