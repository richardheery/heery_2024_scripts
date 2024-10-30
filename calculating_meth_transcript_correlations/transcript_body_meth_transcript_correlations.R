# Calculate methylation-transcription correlations for CPGEA normal, CPGEA tumour and MCRPC samples

# Load required packages
library(methodical)
library(dplyr)

# Get TSS Granges
tss_gr = readRDS("../auxillary_data/pc_transcripts_tss_gr_cage_supported.rds")
transcripts_gr = readRDS("../auxillary_data/pc_transcripts_gr_cage_supported.rds")

# Expand transcripts_gr
transcripts_gr = methodical:::expand_granges(transcripts_gr, 5000, 5000)

# Load CPGEA methylation RSE and transcript counts
cpgea_meth_rse = HDF5Array::loadHDF5SummarizedExperiment("../auxillary_data/methylation_data/cpgea_meth_rse")
cpgea_kallisto_deseq2_counts = data.frame(data.table::fread("../auxillary_data/rnaseq_data/cpgea_kallisto_counts_deseq2_normalized.tsv.gz"), row.names = 1)

# Get CPGEA normal and tumour samples
normal_samples = grep("N", intersect(names(cpgea_kallisto_deseq2_counts), colnames(cpgea_meth_rse)), value = T)
tumour_samples = grep("T", intersect(names(cpgea_kallisto_deseq2_counts), colnames(cpgea_meth_rse)), value = T)

# Create a bpparm object
bpparam = BiocParallel::MulticoreParam(workers = 10) 

# Calculate methylation-transcription correlations for CPGEA normal samples. Took 2.5 hours with 10 cores.  
system.time({transcript_meth_cors_cpgea_normal_samples_5kb = calculateMethSiteTranscriptCors(meth_rse = cpgea_meth_rse, 
  transcript_expression_table = cpgea_kallisto_deseq2_counts, samples_subset = normal_samples, tss_gr = tss_gr, tss_associated_gr = transcripts_gr, 
  cor_method = "spearman", BPPARAM = bpparam, add_distance_to_region = T)})
saveRDS(transcript_meth_cors_cpgea_normal_samples_5kb, "cpgea_normal_whole_gene_body_correlations.rds")

# Calculate methylation-transcription correlations for CPGEA tumour samples. Took 2.5 hours with 10 cores.  
system.time({transcript_meth_cors_cpgea_tumour_samples_5kb = calculateMethSiteTranscriptCors(meth_rse = cpgea_meth_rse, 
  transcript_expression_table = cpgea_kallisto_deseq2_counts, samples_subset = tumour_samples, tss_gr = tss_gr, tss_associated_gr = transcripts_gr, 
  cor_method = "spearman", BPPARAM = bpparam, add_distance_to_region = T)})
saveRDS(transcript_meth_cors_cpgea_tumour_samples_5kb, "cpgea_tumour_whole_gene_body_correlations.rds")

### Calculate correlations for MCRPC
 
# Load MCRPC methylation RSE and transcript counts
mcrpc_meth_rse = HDF5Array::loadHDF5SummarizedExperiment("../auxillary_data/methylation_data/mcrpc_wgbs_hg38")
mcrpc_kallisto_deseq2_counts = data.frame(data.table::fread("../auxillary_data/rnaseq_data/mcrpc_kallisto_counts_deseq2_normalized.tsv.gz"), row.names = 1)

# Get mcrpc normal and tumour samples
common_mcrpc_samples = intersect(names(mcrpc_kallisto_deseq2_counts), colnames(mcrpc_meth_rse))

# Calculate methylation-transcription correlations for MCRPC samples. Took 2.5 hours with 10 cores.  
system.time({transcript_meth_cors_mcrpc_samples_5kb = calculateMethSiteTranscriptCors(meth_rse = mcrpc_meth_rse, 
  transcript_expression_table = mcrpc_kallisto_deseq2_counts, samples_subset = common_mcrpc_samples, tss_gr = tss_gr, tss_associated_gr = transcripts_gr, 
  cor_method = "spearman", BPPARAM = bpparam, add_distance_to_region = T)})
saveRDS(transcript_meth_cors_mcrpc_samples_5kb, "mcrpc_whole_gene_body_correlations.rds")