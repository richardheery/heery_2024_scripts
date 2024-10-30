### Identify TMRs in CPGEA normal, CPGEA tumour and MCRPC datasets

# Load required packages
library(methodical)
library(dplyr)
library(ggplot2)

# Create a BPPARAM object
bpparam = BiocParallel::MulticoreParam(20)

# Load methylation-correlation results for CPGEA normal. Took 3 minutes. 
system.time({cpgea_normal_cor_results = readRDS("../meth_transcript_correlations/cpgea_normal_whole_gene_body_correlations.rds")})

# Find TMRs for CPGEA normal samples. Took 13.5 minutes. 
system.time({cpgea_normal_tmrs = findTMRs(correlation_list = cpgea_normal_cor_results, 
  p_adjust_method = "fdr", p_value_threshold = 0.05, BPPARAM = bpparam)})
rm(cpgea_normal_cor_results); gc()

# Load methylation-correlation results for CPGEA tumour. Took 3 minutes. 
system.time({cpgea_tumour_cor_results = readRDS("../meth_transcript_correlations/cpgea_tumour_whole_gene_body_correlations.rds")})

# Find TMRs for CPGEA tumour samples. Took 13.5 minutes. 
system.time({cpgea_tumour_tmrs = findTMRs(correlation_list = cpgea_tumour_cor_results, 
  p_adjust_method = "fdr", p_value_threshold = 0.05, BPPARAM = bpparam)})
rm(cpgea_tumour_cor_results); gc()

# Load methylation-correlation results for MCRPC. Took 3 minutes. 
system.time({mcrpc_cor_results = readRDS("../meth_transcript_correlations/mcrpc_whole_gene_body_correlations.rds")})

# Find TMRs for MCRPC samples. Took 13.5 minutes. 
system.time({mcrpc_tmrs = findTMRs(correlation_list = mcrpc_cor_results, 
  p_adjust_method = "fdr", p_value_threshold = 0.05, BPPARAM = bpparam)})
rm(mcrpc_cor_results); gc()

# Get repeat ranges for hg38
repeat_ranges = readRDS("~/genomes/repetitive_sequences/repeatmasker/repeatmasker_granges_ucsc.rds")

# Save all TMRs without repeats
cpgea_normal_tmrs_no_repeats = subsetByOverlaps(cpgea_normal_tmrs, repeat_ranges, invert = T)
cpgea_tumour_tmrs_no_repeats = subsetByOverlaps(cpgea_tumour_tmrs, repeat_ranges, invert = T)
mcrpc_tmrs_no_repeats = subsetByOverlaps(mcrpc_tmrs, repeat_ranges, invert = T)

# Create a list of TMRs which do not overlap repeats and save
tmr_list_no_repeats = list(
  cpgea_normal_tmrs_negative = cpgea_normal_tmrs_no_repeats[cpgea_normal_tmrs_no_repeats$direction == "Negative"],
  cpgea_normal_tmrs_positive = cpgea_normal_tmrs_no_repeats[cpgea_normal_tmrs_no_repeats$direction == "Positive"],
  cpgea_tumour_tmrs_negative = cpgea_tumour_tmrs_no_repeats[cpgea_tumour_tmrs_no_repeats$direction == "Negative"],
  cpgea_tumour_tmrs_positive = cpgea_tumour_tmrs_no_repeats[cpgea_tumour_tmrs_no_repeats$direction == "Positive"],
  mcrpc_tmrs_negative = mcrpc_tmrs_no_repeats[mcrpc_tmrs_no_repeats$direction == "Negative"],
  mcrpc_tmrs_positive = mcrpc_tmrs_no_repeats[mcrpc_tmrs_no_repeats$direction == "Positive"]
)
saveRDS(tmr_list_no_repeats, "tmr_list.rds")

# Save CAGE-supported TMRs without repeats
saveRDS(cpgea_normal_tmrs_no_repeats, "gene_body_tmrs/cpgea_normal_tmrs.rds")
saveRDS(cpgea_tumour_tmrs_no_repeats, "gene_body_tmrs/cpgea_tumour_tmrs.rds")
saveRDS(mcrpc_tmrs_no_repeats, "gene_body_tmrs/mcrpc_tmrs.rds")