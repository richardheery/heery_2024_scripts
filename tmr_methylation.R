# Calculate correlation between transcription factor expression and TMR methylation 

# Load required packages
library(methodical)

# Load TMRs
cpgea_normal_tmrs = readRDS("../final_tmr_granges/cpgea_normal_tmrs_5kb.rds")
cpgea_tumour_tmrs = readRDS("../final_tmr_granges/cpgea_tumour_tmrs_5kb.rds")
mcrpc_tmrs = readRDS("../final_tmr_granges/mcrpc_tmrs_5kb.rds")

# Load CPGEA methylation RSE and MCRPC RSE
cpgea_meth_rse = HDF5Array::loadHDF5SummarizedExperiment("~/wgbs/cpgea/wgbs/cpgea_meth_rse_hdf5")
mcrpc_methrix = HDF5Array::loadHDF5SummarizedExperiment("~/wgbs/mcrpc/wgbs/mcrpc_methrix_h5/")
mcrpc_meth_rse = methrix_to_rse(mcrpc_methrix, assays = c("beta"))

# Sort the rowRanges of mcrpc_meth_rse
mcrpc_meth_rse = sort(mcrpc_meth_rse)

# Remove the colData of cpgea_meth_rse and mcrpc_meth_rse so that they can be combined
colData(cpgea_meth_rse) = NULL
colData(mcrpc_meth_rse) = NULL

# Combine cpgea_meth_rse and mcrpc_meth_rse. 
system.time({combined_meth_rse = cbind(cpgea_meth_rse, mcrpc_meth_rse)})

# Get methylation of CPGEA normal TMRs from combined_meth_rse. Took 11 minutes. 
system.time({cpgea_normal_tmrs_meth = summarize_region_methylation(meth_rse = combined_meth_rse, assay_number = 1, 
  genomic_regions = cpgea_normal_tmrs, genomic_regions_names = cpgea_normal_tmrs$tmr_name)})
helpR::fwrite(cpgea_normal_tmrs_meth, "tmr_methylation_tables/cpgea_normal_tmrs_meth_all_samples.tsv.gz")
rm(cpgea_normal_tmrs_meth); gc()

# Get methylation of CPGEA tumour TMRs from combined_meth_rse. Took 11 minutes. 
system.time({cpgea_tumour_tmrs_meth = summarize_region_methylation(meth_rse = combined_meth_rse, assay_number = 1, 
  genomic_regions = cpgea_tumour_tmrs, genomic_regions_names = cpgea_tumour_tmrs$tmr_name)})
helpR::fwrite(cpgea_tumour_tmrs_meth, "tmr_methylation_tables/cpgea_tumour_tmrs_meth_all_samples.tsv.gz")
rm(cpgea_tumour_tmrs_meth); gc()

# Get methylation of MCRPC TMRs from combined_meth_rse. Took 11 minutes. 
system.time({mcrpc_tmrs_meth = summarize_region_methylation(meth_rse = combined_meth_rse, assay_number = 1, 
  genomic_regions = mcrpc_tmrs, genomic_regions_names = mcrpc_tmrs$tmr_name)})
helpR::fwrite(mcrpc_tmrs_meth, "tmr_methylation_tables/mcrpc_tmrs_meth_all_samples.tsv.gz")
