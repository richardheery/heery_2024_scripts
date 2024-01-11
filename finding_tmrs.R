# Find TMRs for CPGEA normal, CPGEA tumour and MCRPC samples

# Make and register a cluster
library(methodical)
library(doParallel)
cl = makeCluster(40)
registerDoParallel(cl, cores = 40)

# Load meth-transcript correlations within 5KB for normal samples
transcript_meth_cors_cpgea_normal_samples = readRDS("../meth_transcript_correlations/meth_transcript_cor_lists/transcript_meth_cors_cpgea_normal_samples_5kb.rds")
system.time({transcript_meth_cors_cpgea_normal_samples = lapply(transcript_meth_cors_cpgea_normal_samples, function(x) dplyr::rename(x, "meth_site" = cpg_name))})

# Run Methodical on normal sample correlations without filtering for CpGs. Took 7 minutes.
system.time({cpgea_normal_tmrs_no_filter = foreach(transcript = transcript_meth_cors_cpgea_normal_samples, .packages = "methodical") %dopar% {
  methodical::find_tmrs(correlation_df = transcript, offset_length = 10, smoothing_factor = 0.75, p_value_threshold = 0.005, min_gapwidth = 150, min_meth_sites = 1)
}})

# Combine all TMRs into a single GRanges and save
system.time({cpgea_normal_tmrs_no_filter = unlist(GRangesList(cpgea_normal_tmrs_no_filter[lengths(cpgea_normal_tmrs_no_filter) > 0]))})
saveRDS(cpgea_normal_tmrs_no_filter, "tmr_granges/cpgea_normal_tmrs_no_filter_5kb.rds")

# Run Methodical on normal sample correlations. Took 7 minutes. 
system.time({cpgea_normal_tmrs = foreach(transcript = transcript_meth_cors_cpgea_normal_samples, .packages = "methodical") %dopar% {
  methodical::find_tmrs(correlation_df = transcript, offset_length = 10, smoothing_factor = 0.75, p_value_threshold = 0.005, min_gapwidth = 150, min_meth_sites = 5)
}})

# Combine all TMRs into a single GRanges and save
system.time({cpgea_normal_tmrs = unlist(GRangesList(cpgea_normal_tmrs[lengths(cpgea_normal_tmrs) > 0]))})
saveRDS(cpgea_normal_tmrs, "tmr_granges/cpgea_normal_tmrs_5kb.rds")

# Remove transcript_meth_cors_cpgea_normal_samples
rm(transcript_meth_cors_cpgea_normal_samples); gc()

# Load meth-transcript correlations within 5KB for tumour samples
transcript_meth_cors_cpgea_tumour_samples = readRDS("../meth_transcript_correlations/meth_transcript_cor_lists/transcript_meth_cors_cpgea_tumour_samples_5kb.rds")
system.time({transcript_meth_cors_cpgea_tumour_samples = lapply(transcript_meth_cors_cpgea_tumour_samples, function(x) dplyr::rename(x, "meth_site" = cpg_name))})

# Run Methodical on tumour sample correlations. Took 7 minutes. 
system.time({cpgea_tumour_tmrs = foreach(transcript = transcript_meth_cors_cpgea_tumour_samples, .packages = "methodical") %dopar% {
  find_tmrs(correlation_df = transcript, offset_length = 10, smoothing_factor = 0.75, p_value_threshold = 0.005, min_gapwidth = 150, min_meth_sites = 5)
}})

# Combine all TMRs into a single GRanges and save
system.time({cpgea_tumour_tmrs = unlist(GRangesList(cpgea_tumour_tmrs[lengths(cpgea_tumour_tmrs) > 0]))})
saveRDS(cpgea_tumour_tmrs, "tmr_granges/cpgea_tumour_tmrs_5kb.rds")

# Remove transcript_meth_cors_cpgea_normal_samples
rm(transcript_meth_cors_cpgea_tumour_samples); gc()

# Load meth-transcript correlations within 5KB for metastasis samples
transcript_meth_cors_mcrpc_samples = readRDS("../meth_transcript_correlations/meth_transcript_cor_lists/transcript_meth_cors_mcrpc_samples_5kb.rds")
system.time({transcript_meth_cors_mcrpc_samples = lapply(transcript_meth_cors_mcrpc_samples, function(x) dplyr::rename(x, "meth_site" = cpg_name))})

# Run Methodical on MCRPC sample correlations. Took 7 minutes. 
system.time({mcrpc_tmrs = foreach(transcript = transcript_meth_cors_mcrpc_samples, .packages = "methodical") %dopar% {
  find_tmrs(correlation_df = transcript, offset_length = 10, smoothing_factor = 0.75, p_value_threshold = 0.005, min_gapwidth = 150, min_meth_sites = 5)
}})

# Combine all TMRs into a single GRanges and save
system.time({mcrpc_tmrs = unlist(GRangesList(mcrpc_tmrs[lengths(mcrpc_tmrs) > 0]))})
saveRDS(mcrpc_tmrs, "tmr_granges/mcrpc_tmrs_5kb.rds")

## Load meth-transcript correlations within 5KB for rhabdoid
transcript_meth_cors_rhabdoid_samples = readRDS("../meth_transcript_correlations/meth_transcript_cor_lists/transcript_meth_cors_rhabdoid_samples_5kb.rds")

# Run Methodical on rhabdoid sample correlations. Took 6 minutes. 
system.time({rhabdoid_tmrs = foreach(transcript = transcript_meth_cors_rhabdoid_samples, .packages = "methodical") %dopar% {
  find_tmrs(correlation_df = transcript, offset_length = 10, smoothing_factor = 0.75, p_value_threshold = 0.005, min_gapwidth = 150, min_meth_sites = 5)
}})

# Combine all TMRs into a single GRanges and save
system.time({rhabdoid_tmrs = unlist(GRangesList(rhabdoid_tmrs[lengths(rhabdoid_tmrs) > 0]))})
saveRDS(rhabdoid_tmrs, "tmr_granges/rhabdoid_tmrs_5kb.rds")

### Find TMRs in +/- 50 KB and +/- 500 KB

# Find 50 KB TMRs for CPGEA normal samples. Took 7 minutes with 35 cores. 
transcript_meth_cors_cpgea_normal_50kb_files = list.files("../meth_transcript_correlations//meth_transcript_cor_directories/transcript_meth_cors_cpgea_normal_samples_50kb_files/", full.names = T) 
system.time({cpgea_normal_tmrs_50kb = foreach(transcript = transcript_meth_cors_cpgea_normal_50kb_files, .packages = "methodical") %dopar% {
  find_tmrs(correlation_df = transcript, offset_length = 10, smoothing_factor = 0.75, p_value_threshold = 0.005)}})

# Combine all TMRs into a single GRanges and save. Took 2 minutes. 
system.time({cpgea_normal_tmrs_50kb = unlist(GRangesList(cpgea_normal_tmrs_50kb[lengths(cpgea_normal_tmrs_50kb) > 0]))})
saveRDS(cpgea_normal_tmrs_50kb, "tmr_granges/cpgea_normal_tmrs_50kb.rds")

# Find 50 KB TMRs for CPGEA tumour samples. Took 7 minutes with 35 cores. 
transcript_meth_cors_cpgea_tumour_50kb_files = list.files("../meth_transcript_correlations//meth_transcript_cor_directories/transcript_meth_cors_cpgea_tumour_samples_50kb_files/", full.names = T) 
system.time({cpgea_tumour_tmrs_50kb = foreach(transcript = transcript_meth_cors_cpgea_tumour_50kb_files, .packages = "methodical") %dopar% {
  find_tmrs(correlation_df = transcript, offset_length = 10, smoothing_factor = 0.75, p_value_threshold = 0.005)}})

# Combine all TMRs into a single GRanges and save. Took 3 minutes. 
system.time({cpgea_tumour_tmrs_50kb = unlist(GRangesList(cpgea_tumour_tmrs_50kb[lengths(cpgea_tumour_tmrs_50kb) > 0]))})
saveRDS(cpgea_tumour_tmrs_50kb, "tmr_granges/cpgea_tumour_tmrs_50kb.rds")

# Find 50 KB TMRs for mcrpc samples. Took 7 minutes with 35 cores. 
transcript_meth_cors_mcrpc_50kb_files = list.files("../meth_transcript_correlations//meth_transcript_cor_directories/transcript_meth_cors_mcrpc_samples_50kb_files/", full.names = T) 
system.time({mcrpc_tmrs_50kb = foreach(transcript = transcript_meth_cors_mcrpc_50kb_files, .packages = "methodical") %dopar% {
  find_tmrs(correlation_df = transcript, offset_length = 10, smoothing_factor = 0.75, p_value_threshold = 0.005)}})

# Combine all TMRs into a single GRanges and save. Took 2 minutes. 
system.time({mcrpc_tmrs_50kb = unlist(GRangesList(mcrpc_tmrs_50kb[lengths(mcrpc_tmrs_50kb) > 0]))})
saveRDS(mcrpc_tmrs_50kb, "tmr_granges/mcrpc_tmrs_50kb.rds")

# Find 500 KB TMRs for CPGEA normal samples. Took 37 minutes with 35 cores. 
transcript_meth_cors_cpgea_normal_500kb_files = list.files("../meth_transcript_correlations//meth_transcript_cor_directories/transcript_meth_cors_cpgea_normal_samples_500kb_files/", full.names = T) 
system.time({cpgea_normal_tmrs_500kb = foreach(transcript = transcript_meth_cors_cpgea_normal_500kb_files, .packages = "methodical") %dopar% {
  find_tmrs(correlation_df = transcript, offset_length = 10, smoothing_factor = 0.75, p_value_threshold = 0.005)}})

# Combine all TMRs into a single GRanges and save. Took 5 minutes. 
system.time({cpgea_normal_tmrs_500kb = unlist(GRangesList(cpgea_normal_tmrs_500kb[lengths(cpgea_normal_tmrs_500kb) > 0]))})
saveRDS(cpgea_normal_tmrs_500kb, "tmr_granges/cpgea_normal_tmrs_500kb.rds")

# Find 500 KB TMRs for CPGEA tumour samples. Took 36 minutes with 35 cores. 
transcript_meth_cors_cpgea_tumour_500kb_files = list.files("../meth_transcript_correlations//meth_transcript_cor_directories/transcript_meth_cors_cpgea_tumour_samples_500kb_files/", full.names = T) 
system.time({cpgea_tumour_tmrs_500kb = foreach(transcript = transcript_meth_cors_cpgea_tumour_500kb_files, .packages = "methodical") %dopar% {
  find_tmrs(correlation_df = transcript, offset_length = 10, smoothing_factor = 0.75, p_value_threshold = 0.005)}})

# Combine all TMRs into a single GRanges and save. Took 5 minutes. 
system.time({cpgea_tumour_tmrs_500kb = unlist(GRangesList(cpgea_tumour_tmrs_500kb[lengths(cpgea_tumour_tmrs_500kb) > 0]))})
saveRDS(cpgea_tumour_tmrs_500kb, "tmr_granges/cpgea_tumour_tmrs_500kb.rds")

# Find 500 KB TMRs for MCRPC samples. Took 34 minutes with 35 cores. 
transcript_meth_cors_mcrpc_500kb_files = list.files("../meth_transcript_correlations//meth_transcript_cor_directories/transcript_meth_cors_mcrpc_samples_500kb_files/", full.names = T) 
system.time({mcrpc_tmrs_500kb = foreach(transcript = transcript_meth_cors_mcrpc_500kb_files, .packages = "methodical") %dopar% {
  find_tmrs(correlation_df = transcript, offset_length = 10, smoothing_factor = 0.75, p_value_threshold = 0.005)}})

# Combine all TMRs into a single GRanges and save. Took 5 minutes. 
system.time({mcrpc_tmrs_500kb = unlist(GRangesList(mcrpc_tmrs_500kb[lengths(mcrpc_tmrs_500kb) > 0]))})
saveRDS(mcrpc_tmrs_500kb, "tmr_granges/mcrpc_tmrs_500kb.rds")

### Filter TMRs for those not overlapping repeats

# Get paths to all confident (those with at least 5 CpGs) TMR GRanges
tmr_granges_files = list.files("tmr_granges", full.names = T, pattern = "tmrs_5")

# Get repeat ranges for hg38
repeat_ranges = readRDS("~/genomes/repetitive_sequences/repeatmasker/repeatmasker_granges_ucsc.rds")

# Read in each file, filter filter for TMRs not overlapping repeats and save to final_tmr_granges
foreach(tmr_file = tmr_granges_files) %do% {
  tmrs = readRDS(tmr_file)
  tmrs = IRanges::subsetByOverlaps(tmrs, repeat_ranges, invert = T)
  saveRDS(tmrs, paste0("final_tmr_granges/", basename(tmr_file)))
}

# Get paths to all Final (those not overlapping repeats) TMR GRanges
final_tmr_granges_files = list.files("final_tmr_granges", full.names = T)

# Get paths to 5KB TMR files for each dataset
final_tmr_5kb_files = list.files("final_tmr_granges", full.names = T, pattern = "tmrs_5kb")
names(final_tmr_5kb_files) = gsub("_5kb", "", basename(tools::file_path_sans_ext(final_tmr_5kb_files)))

# Create a list with 5KB final_tmrs for each dataset
final_tmr_5kb_list = lapply(final_tmr_5kb_files, function(x) readRDS(x))

# Split groups into negative and positive
final_tmr_5kb_list = unlist(lapply(final_tmr_5kb_list, function(x) split(x, x$direction)))
names(final_tmr_5kb_list) = tolower(gsub("\\.", "_", names(final_tmr_5kb_list)))
saveRDS(final_tmr_5kb_list, "final_tmr_granges/final_tmr_5kb_list.rds")

### Test enrichment of TMRs for transcripts associated with specific tags

# Load all TMRs
final_tmr_5kb_list = readRDS("final_tmr_granges/final_tmr_5kb_list.rds")

# Create BED files for final_tmrs
lapply(names(final_tmr_5kb_list), function(x)
  rtracklayer::export.bed(object =  setNames(final_tmr_5kb_list[[x]], final_tmr_5kb_list[[x]]$tmr_name), paste0("final_tmr_granges/bed_files/", x, ".bed")))

# Get the IDs for protein-coding transcripts
pcg_transcripts = readRDS("~/genomes/gencode/gencode_granges/pcg_transcript_ranges_gencode_v38.rds")$transcript_id

# Get the GFF3 Gencode annotation
gencode_annotation = rtracklayer::import.gff3("~/genomes/gencode/gencode_downloads/gencode_38/gencode.v38.annotation.gff3.gz")

# Filter for transcript annotation
gencode_annotation = gencode_annotation[gencode_annotation$type == "transcript"]

# Get tags associated with each transcript
gencode_transcript_tags = as.list(setNames(gencode_annotation$tag, gsub("\\.[0-9]*", "", gencode_annotation$ID)))
gencode_transcript_tags = gencode_transcript_tags[pcg_transcripts]

# Find the proportion of transcripts with each tag
gencode_transcript_tags_counts = table(unlist(gencode_transcript_tags))/length(gencode_transcript_tags)

# Get unique tags 
tags = sort(unique(names(gencode_transcript_tags_counts)))
names(tags) = tags

# For each TMR group, check which transcripts are associated with TMRs
system.time({transcripts_with_tmrs = data.frame(lapply(final_tmr_5kb_list, function(x) 
    setNames(pcg_transcripts, pcg_transcripts) %in% x$transcript_id))})

# For each tag, check which transcripts are associated with it
system.time({transcripts_with_tags = data.frame(lapply(tags, function(x) 
    sapply(gencode_transcript_tags, function(y) x %in% y)))})

# Perform Fisher tests testing if TMR-associated transcripts tend to be also associated with certain tags
fisher_test_results = bind_rows(lapply(transcripts_with_tmrs, function(x)
    bind_rows(lapply(transcripts_with_tags, function(y)
        broom::tidy(fisher.test(x, y))), .id = "tag")), .id = "tmr_group")

# Add q-value to the result
fisher_test_results$q_value = p.adjust(fisher_test_results$p.value, method = "fdr")

# Remove X from start of 3' or 5' related tags
fisher_test_results$tag = gsub("^X", "", fisher_test_results$tag)

# Get significantly enriched results
# No significant results for cpgea_normal_tmrs_positive
# MANE select, Ensembl Canonical, APPRIS principal 1, Gencode basic and CCDS enriched for all TMR groups except normal positive. 
# MANE select most enriched for each of the 5 groups, with Ensembl canonical 2nd. Both enriched 2-4 times depending on TMR group. 
# All negative TMRs enriched for CAGE supported TSS
fisher_test_results_significant = arrange(filter(fisher_test_results, q_value < 0.05, estimate > 1), q_value)
