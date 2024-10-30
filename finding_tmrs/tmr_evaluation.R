# Evaluate the TMR correlations

# Load required packages
library(methodical)
library(dplyr)
library(patchwork)
source("../auxillary_scripts/plotting_functions.R")

# Load CPGEA methylation RSE and MCRPC RSE
cpgea_meth_rse = HDF5Array::loadHDF5SummarizedExperiment("../auxillary_data/methylation_data/cpgea_wgbs_hg38")
mcrpc_meth_rse =  HDF5Array::loadHDF5SummarizedExperiment("../auxillary_data/methylation_data/mcrpc_wgbs_hg38")

# Get transcript counts for CPGEA and MCRPC
cpgea_kallisto_deseq2_counts = data.frame(data.table::fread("../auxillary_data/rnaseq_data/cpgea_kallisto_counts_deseq2_normalized.tsv.gz"), row.names = 1)
mcrpc_kallisto_deseq2_counts = data.frame(data.table::fread("../auxillary_data/rnaseq_data/mcrpc_kallisto_counts_deseq2_normalized.tsv.gz"), row.names = 1)

# Get CPGEA normal and tumour samples and MCRPC samples with both WGBS and RNA-seq data
common_cpgea_normal_samples = grep("N", intersect(names(cpgea_kallisto_deseq2_counts), colnames(cpgea_meth_rse)), value = T)
common_cpgea_tumour_samples = grep("T", intersect(names(cpgea_kallisto_deseq2_counts), colnames(cpgea_meth_rse)), value = T)
common_mcrpc_samples = intersect(names(mcrpc_kallisto_deseq2_counts), colnames(mcrpc_meth_rse))

# Get TMR GRanges
cpgea_normal_tmrs = readRDS("cpgea_normal_tmrs.rds")
cpgea_tumour_tmrs = readRDS("cpgea_tumour_tmrs.rds")
mcrpc_tmrs = readRDS("mcrpc_tmrs.rds")

### Evaluate TMR correlations in normal prostate samples

# Took 1 minutes with 10 cores
system.time({cpgea_normal_tmrs_correlations_normal_samples = calculateRegionMethylationTranscriptCors(
  genomic_regions = cpgea_normal_tmrs, genomic_region_names = cpgea_normal_tmrs$tmr_name, samples_subset = common_cpgea_normal_samples,
  genomic_region_transcripts = cpgea_normal_tmrs$ID, meth_rse = cpgea_meth_rse, transcript_expression_table = cpgea_kallisto_deseq2_counts, 
  cor_method = "s", ncores = 10)})
data.table::fwrite(cpgea_normal_tmrs_correlations_normal_samples, "tmr_evaluation_tables/cpgea_normal_tmrs_correlations_normal_samples.tsv.gz")

# Took 1 minutes with 10 cores
system.time({cpgea_tumour_tmrs_correlations_normal_samples = calculateRegionMethylationTranscriptCors(
  genomic_regions = cpgea_tumour_tmrs, genomic_region_names = cpgea_tumour_tmrs$tmr_name, samples_subset = common_cpgea_normal_samples,
  genomic_region_transcripts = cpgea_tumour_tmrs$ID, meth_rse = cpgea_meth_rse, transcript_expression_table = cpgea_kallisto_deseq2_counts, 
  cor_method = "s", ncores = 10)})
data.table::fwrite(cpgea_tumour_tmrs_correlations_normal_samples, "tmr_evaluation_tables/cpgea_tumour_tmrs_correlations_normal_samples.tsv.gz")

# Took 1 minutes with 10 cores
system.time({mcrpc_tmrs_correlations_normal_samples = calculateRegionMethylationTranscriptCors(
  genomic_regions = mcrpc_tmrs, genomic_region_names = mcrpc_tmrs$tmr_name, samples_subset = common_cpgea_normal_samples,
  genomic_region_transcripts = mcrpc_tmrs$ID, meth_rse = cpgea_meth_rse, transcript_expression_table = cpgea_kallisto_deseq2_counts, 
  cor_method = "s", ncores = 10)})
data.table::fwrite(mcrpc_tmrs_correlations_normal_samples, "tmr_evaluation_tables/mcrpc_tmrs_correlations_normal_samples.tsv.gz")

### Evaluate TMR correlations in prostate tumour samples

# Took 1 minutes with 10 cores
system.time({cpgea_normal_tmrs_correlations_tumour_samples = calculateRegionMethylationTranscriptCors(
  genomic_regions = cpgea_normal_tmrs, genomic_region_names = cpgea_normal_tmrs$tmr_name, samples_subset = common_cpgea_tumour_samples,
  genomic_region_transcripts = cpgea_normal_tmrs$ID, meth_rse = cpgea_meth_rse, transcript_expression_table = cpgea_kallisto_deseq2_counts, 
  cor_method = "s", ncores = 10)})
data.table::fwrite(cpgea_normal_tmrs_correlations_tumour_samples, "tmr_evaluation_tables/cpgea_normal_tmrs_correlations_tumour_samples.tsv.gz")

# Took 1 minutes with 10 cores
system.time({cpgea_tumour_tmrs_correlations_tumour_samples = calculateRegionMethylationTranscriptCors(
  genomic_regions = cpgea_tumour_tmrs, genomic_region_names = cpgea_tumour_tmrs$tmr_name, samples_subset = common_cpgea_tumour_samples,
  genomic_region_transcripts = cpgea_tumour_tmrs$ID, meth_rse = cpgea_meth_rse, transcript_expression_table = cpgea_kallisto_deseq2_counts, 
  cor_method = "s", ncores = 10)})
data.table::fwrite(cpgea_tumour_tmrs_correlations_tumour_samples, "tmr_evaluation_tables/cpgea_tumour_tmrs_correlations_tumour_samples.tsv.gz")

# Took 1 minutes with 10 cores
system.time({mcrpc_tmrs_correlations_tumour_samples = calculateRegionMethylationTranscriptCors(
  genomic_regions = mcrpc_tmrs, genomic_region_names = mcrpc_tmrs$tmr_name, samples_subset = common_cpgea_tumour_samples,
  genomic_region_transcripts = mcrpc_tmrs$ID, meth_rse = cpgea_meth_rse, transcript_expression_table = cpgea_kallisto_deseq2_counts, 
  cor_method = "s", ncores = 10)})
data.table::fwrite(mcrpc_tmrs_correlations_tumour_samples, "tmr_evaluation_tables/mcrpc_tmrs_correlations_tumour_samples.tsv.gz")

### Evaluate TMR correlations in metastases samples

# Took 1 minutes with 10 cores
system.time({cpgea_normal_tmrs_correlations_metastases_samples = calculateRegionMethylationTranscriptCors(
  genomic_regions = cpgea_normal_tmrs, genomic_region_names = cpgea_normal_tmrs$tmr_name, samples_subset = common_mcrpc_samples,
  genomic_region_transcripts = cpgea_normal_tmrs$ID, meth_rse = mcrpc_meth_rse, transcript_expression_table = mcrpc_kallisto_deseq2_counts, 
  cor_method = "s", ncores = 10)})
data.table::fwrite(cpgea_normal_tmrs_correlations_metastases_samples, "tmr_evaluation_tables/cpgea_normal_tmrs_correlations_metastases_samples.tsv.gz")

# Took 1 minutes with 10 cores
system.time({cpgea_tumour_tmrs_correlations_metastases_samples = calculateRegionMethylationTranscriptCors(
  genomic_regions = cpgea_tumour_tmrs, genomic_region_names = cpgea_tumour_tmrs$tmr_name, samples_subset = common_mcrpc_samples,
  genomic_region_transcripts = cpgea_tumour_tmrs$ID, meth_rse = mcrpc_meth_rse, transcript_expression_table = mcrpc_kallisto_deseq2_counts, 
  cor_method = "s", ncores = 10)})
data.table::fwrite(cpgea_tumour_tmrs_correlations_metastases_samples, "tmr_evaluation_tables/cpgea_tumour_tmrs_correlations_metastases_samples.tsv.gz")

# Took 1 minutes with 10 cores
system.time({mcrpc_tmrs_correlations_metastases_samples = calculateRegionMethylationTranscriptCors(
  genomic_regions = mcrpc_tmrs, genomic_region_names = mcrpc_tmrs$tmr_name, samples_subset = common_mcrpc_samples,
  genomic_region_transcripts = mcrpc_tmrs$ID, meth_rse = mcrpc_meth_rse, transcript_expression_table = mcrpc_kallisto_deseq2_counts, 
  cor_method = "s", ncores = 10)})
data.table::fwrite(mcrpc_tmrs_correlations_metastases_samples, "tmr_evaluation_tables/mcrpc_tmrs_correlations_metastases_samples.tsv.gz")

### Create plots of cluster correlations

# Get paths to all cluster evaluation tables
cluster_evaluation_tables_paths = list.files("tmr_evaluation_tables/", full.names = T)
names(cluster_evaluation_tables_paths) = gsub(".tsv.gz", "", basename(cluster_evaluation_tables_paths))

# Read in all tables as a list
cluster_evaluation_tables = lapply(cluster_evaluation_tables_paths, function(x) data.table::fread(x))

# Add direction to TMRs
for(evaluation in names(cluster_evaluation_tables)){

  cluster_group = get(gsub("_correlations.*", "", evaluation))
  cluster_evaluation_tables[[evaluation]]$direction = cluster_group$direction[match(
    cluster_evaluation_tables[[evaluation]]$genomic_region_name, cluster_group$tmr_name)]

}

# Combine all tables 
cluster_evaluation_tables = bind_rows(cluster_evaluation_tables, .id = "evaluation")

# Update q-values
cluster_evaluation_tables$q_val = p.adjust(cluster_evaluation_tables$p_val, method = "fdr")

# Add TMR group and sample group to table
cluster_evaluation_tables$tmr_group = gsub("_correlations.*", "", cluster_evaluation_tables$evaluation)
cluster_evaluation_tables$samples = gsub(".*_correlations_", "", cluster_evaluation_tables$evaluation)

# Add a column indicating which samples the TMRs are from so that external correlations can be easily identified
cluster_evaluation_tables$tmr_samples = case_when(
  cluster_evaluation_tables$tmr_group == "cpgea_normal_tmrs" ~ "normal_samples",
  cluster_evaluation_tables$tmr_group == "cpgea_tumour_tmrs" ~ "tumour_samples",
  cluster_evaluation_tables$tmr_group == "mcrpc_tmrs" ~ "metastases_samples")

# Calculate proportion of significant external correlations first for all datasets and then just considering the tumour and metastases datasets
prop.table(table(filter(cluster_evaluation_tables, samples != tmr_samples)$q_val < 0.05))
prop.table(table(filter(cluster_evaluation_tables, samples != tmr_samples & samples != "normal_samples" & tmr_samples != "normal_samples")$q_val < 0.05))

# Convert samples and tmr_group to factor and put levels in correct order
cluster_evaluation_tables$samples = factor(cluster_evaluation_tables$samples, levels = c("normal_samples", "tumour_samples", "metastases_samples"))
cluster_evaluation_tables$tmr_group = factor(cluster_evaluation_tables$tmr_group, 
  levels = c("cpgea_normal_tmrs", "cpgea_tumour_tmrs", "mcrpc_tmrs"))

# Create violin plots showing the distribution TMR correlations
correlation_distribution_violins = ggplot(cluster_evaluation_tables, aes(x = samples, y = cor, fill = direction)) +
  geom_violin() + geom_hline(yintercept = 0, linetype = "dashed")
correlation_distribution_violins = customize_ggplot_theme(correlation_distribution_violins, facet = "tmr_group", facet_nrow = 3, fill_colors = colour_list$purple_and_gold_light, 
  x_labels = c("Normal Prostate Samples", "Prostate Tumour Samples", "Prostate Metastasis Samples"), 
  ylab = "Methylation-Transcription Correlation", fill_title = "TMR Direction",
  fill_labels = c("Negative", "Positive"), facet_labels = c("Normal Prostate TMRs", "Prostate Tumour TMRs", "Prostate Metastasis TMRs")) +
  theme(strip.background = element_blank())
correlation_distribution_violins

# Count number of significant correlations by direction, tmr_group and sample group
significant_correlation_df = summarize(group_by(cluster_evaluation_tables, samples, tmr_group, direction), 
  proportion_sig = sum(q_val < 0.05, na.rm = T)/dplyr::n())

# Create barplots showing the number of significant correlations
significant_tmr_cor_barplots = ggplot(significant_correlation_df, aes(x = samples, y = proportion_sig, fill = direction)) +
  geom_col(position = "dodge", color = "black")
significant_tmr_cor_barplots = customize_ggplot_theme(significant_tmr_cor_barplots, facet = "tmr_group", facet_nrow = 3, fill_colors = colour_list$purple_and_gold_light, 
  x_labels = c("Normal Prostate Samples", "Prostate Tumour Samples", "Prostate Metastasis Samples"), 
  ylab = "Proportion of Correlations that are Significant", fill_title = "TMR Direction",
  fill_labels = c("Negative", "Positive"), facet_labels = c("Normal Prostate TMRs", "Prostate Tumour TMRs", "Prostate Metastasis TMRs")) +
  scale_y_continuous(limits = c(0, 1), expand = c(0, 0)) +
  theme(strip.background = element_blank())
significant_tmr_cor_barplots

# Combine correlation_distribution_violins and significant_tmr_cor_barplots
combined_violins_and_barplot = correlation_distribution_violins / significant_tmr_cor_barplots
ggsave(plot = combined_violins_and_barplot, "../figures/supplementary_figure8D_and_E.pdf", width = 27, height = 23.63)

# Count proportion of significant correlations for external TMRs 
cluster_evaluation_tables_external = filter(cluster_evaluation_tables, samples != tmr_samples)
sum(cluster_evaluation_tables_external$q_val < 0.05, na.rm = T)/sum(!is.na(cluster_evaluation_tables_external$q_val))

# Count proportion of significant correlations for external TMRs involving prostate tumours and prostate metastases 
cluster_evaluation_tables_external_cancer = filter(cluster_evaluation_tables_external, samples != "normal_samples" & tmr_samples != "normal_samples")
sum(cluster_evaluation_tables_external_cancer$q_val < 0.05, na.rm = T)/sum(!is.na(cluster_evaluation_tables_external_cancer$q_val))