# Create plots comparing WGBS and Probe correlation results

# Load required packages
library(methodical)
library(SummarizedExperiment)
library(dplyr)
source("../auxillary_scripts/plotting_functions.R")

### Load necessary data

# Get Illumina probes ranges for hg38
illumina_450k_probes_hg38 = readRDS("../auxillary_data/genomic_annotation/infinium_450k_probe_granges_hg38.rds")
illumina_epic_probes_hg38 = readRDS("../auxillary_data/genomic_annotation//epic_probe_gr_hg38.rds")

# Get names of TSS supported by CAGE
cage_supported_tss = readRDS("../auxillary_data/genomic_annotation/cage_supported_transcript_ids.rds")

# Get correlation results for CPGEA normal samples and combine 
cpgea_normal_correlations = readRDS("../calculating_meth_transcript_correlations/cpgea_normal_whole_gene_body_correlations.rds")[cage_supported_tss]
cpgea_normal_correlations = dplyr::bind_rows(cpgea_normal_correlations, .id = "transcript_id")

# Correct p-values
cpgea_normal_correlations$q_val = p.adjust(cpgea_normal_correlations$p_val, method = "fdr")

# Bin CpGs by distance to TSS
cpgea_normal_correlations$bin = plyr::round_any(cpgea_normal_correlations$distance_to_tss, 500)

# Add columns indicating if CpGs are covered by 450K and EPIC probes
cpgea_normal_correlations$probe_450k = cpgea_normal_correlations$meth_site %in% as.character(illumina_450k_probes_hg38)
cpgea_normal_correlations$probe_epic = cpgea_normal_correlations$meth_site %in% as.character(illumina_epic_probes_hg38)

### Create plot of the number of CpG sites in bins around TSS

# Count mean number of CpGs in each bin for each transcript. Took 18 minutes
system.time({cpgs_per_bin = summarise(group_by(cpgea_normal_correlations, bin), 
  count = n()/length(unique(cpgea_normal_correlations$transcript_id)))})
cpgs_per_bin = filter(cpgs_per_bin, abs(bin) < 5000)

# Create plot of mean number of CpGs per bin
cpgs_per_bin_plot = ggplot(cpgs_per_bin, aes(x = bin, y = count)) +
  geom_col(fill = "#BCBDDC")
cpgs_per_bin_plot = customize_ggplot_theme(cpgs_per_bin_plot, ylab = "Mean Number of\nCpG Sites", xlab = "Distance to TSS (bp)") +
  scale_x_continuous(expand = c(0, 0), labels = scales::comma)
cpgs_per_bin_plot
saveRDS(cpgs_per_bin_plot, "cpgs_per_bin_plot.rds")

### Create plot of the number of CpG sites targeted by Illumina arrays in bins around TSS

# Count mean number of CpG covered by 450k probes in each bin for each transcript
probes_450k_per_bin = summarise(group_by(filter(cpgea_normal_correlations, probe_450k, abs(bin) < 5000), bin), 
  count = n()/length(unique(cpgea_normal_correlations$transcript_id)))

# Count mean number of CpG covered by epic probes in each bin for each transcript
probes_epic_per_bin = summarise(group_by(filter(cpgea_normal_correlations, probe_epic, abs(bin) < 5000), bin), 
  count = n()/length(unique(cpgea_normal_correlations$transcript_id)))

# Combine probes_450k _per_bin and probes_epic_per_bin
combined_probes_per_bin = bind_rows(list(probes_450k = probes_450k_per_bin, probes_epic = probes_epic_per_bin), .id = "probe")
combined_probes_per_bin = arrange(combined_probes_per_bin, desc(probe))

# Create plot of mean number of probes per bin
combined_probes_per_bin_plot = ggplot(combined_probes_per_bin, aes(x = bin, y = count, fill = probe)) +
  geom_col(position = "dodge", color = "black")
combined_probes_per_bin_plot = customize_ggplot_theme(combined_probes_per_bin_plot, 
  ylab = "Mean Number of\nIllumina Array Probes", xlab = "Distance to TSS (bp)", 
  fill_colors = colour_list$three_blues[c(1, 2)], fill_labels = c("450K Probes", "EPIC Probes")) +
  scale_x_continuous(expand = c(0, 0), labels = scales::comma) +
  theme(legend.position = c(0.85, 0.85))
combined_probes_per_bin_plot
saveRDS(combined_probes_per_bin_plot, "combined_probes_per_bin_plot.rds")

### Create a plot of the proportion of significant correlations in bins in normal samples

# Bin correlations and count number of significant correlations and all correlations in bin. 
cpgea_normal_all_correlations_binned = data.frame(summarize(group_by(cpgea_normal_correlations, bin), 
  num_sig = sum(q_val < 0.05, na.rm = T), total = sum(!is.na(p_val))))

# Calculate proportion of significant correlations
cpgea_normal_all_correlations_binned$prop_sig = cpgea_normal_all_correlations_binned$num_sig/cpgea_normal_all_correlations_binned$total

# Create barplot of significant correlations per bin and the percentage of significant correlations covered by probes
cpgea_normal_prop_sig_bins_plot = ggplot(filter(cpgea_normal_all_correlations_binned, abs(bin) < 5000), aes(x = bin, y = prop_sig)) +
  geom_col(position = "identity", fill = "#2a5674") 
cpgea_normal_prop_sig_bins_plot = customize_ggplot_theme(cpgea_normal_prop_sig_bins_plot, xlab = "Distance to TSS (bp)", 
  ylab = "Proportion of CpG Sites\nDisplaying Significant Correlations", title = NULL)
cpgea_normal_prop_sig_bins_plot = cpgea_normal_prop_sig_bins_plot + scale_x_continuous(expand = c(0, 0), labels = scales::comma)
cpgea_normal_prop_sig_bins_plot
saveRDS(cpgea_normal_prop_sig_bins_plot, "cpgea_normal_prop_sig_bins_plot.rds")

### Create significance plots for tumour samples

# Get correlation results for CPGEA tumour samples and combine them and add column for bin
cpgea_tumour_correlations = readRDS("../calculating_meth_transcript_correlations/cpgea_tumour_whole_gene_body_correlations.rds")[cage_supported_tss]
cpgea_tumour_correlations = dplyr::bind_rows(cpgea_tumour_correlations, .id = "transcript_id")
cpgea_tumour_correlations$bin = plyr::round_any(cpgea_tumour_correlations$distance_to_tss, 500)

# Correct p-values
cpgea_tumour_correlations$q_val = p.adjust(cpgea_tumour_correlations$p_val, method = "fdr")

# Bin correlations and count number of significant correlations and all correlations in bin. 
cpgea_tumour_all_correlations_binned = data.frame(summarize(group_by(cpgea_tumour_correlations, bin), 
  num_sig = sum(q_val < 0.05, na.rm = T), total = sum(!is.na(p_val))))

# Calculate proportion of significant correlations
cpgea_tumour_all_correlations_binned$prop_sig = cpgea_tumour_all_correlations_binned$num_sig/cpgea_tumour_all_correlations_binned$total

# Create barplot of significant correlations per bin and the percentage of significant correlations covered by probes
cpgea_tumour_prop_sig_bins_plot = ggplot(filter(cpgea_tumour_all_correlations_binned, abs(bin) < 5000), aes(x = bin, y = prop_sig)) +
  geom_col(position = "identity", fill = "#2a5674") 
cpgea_tumour_prop_sig_bins_plot = customize_ggplot_theme(cpgea_tumour_prop_sig_bins_plot, xlab = "Distance to TSS (bp)", 
  ylab = "Proportion of CpG Sites\nDisplaying Significant Correlations", title = NULL)
cpgea_tumour_prop_sig_bins_plot = cpgea_tumour_prop_sig_bins_plot + scale_x_continuous(expand = c(0, 0), labels = scales::comma)
cpgea_tumour_prop_sig_bins_plot
saveRDS(cpgea_tumour_prop_sig_bins_plot, "cpgea_tumour_prop_sig_bins_plot.rds")

### Create plots for MCRPC metastasis samples

# Get correlation results for CPGEA normal samples and combine them and add column for bin
mcrpc_correlations = readRDS("../calculating_meth_transcript_correlations/mcrpc_whole_gene_body_correlations.rds")[cage_supported_tss]
mcrpc_correlations = dplyr::bind_rows(mcrpc_correlations, .id = "transcript_id")
mcrpc_correlations$bin = plyr::round_any(mcrpc_correlations$distance_to_tss, 500)

# Correct p-values
mcrpc_correlations$q_val = p.adjust(mcrpc_correlations$p_val, method = "fdr")

# Bin correlations and count number of significant correlations and all correlations in bin. 
# Significant is defined as a q-value < 0.05
mcrpc_all_correlations_binned = data.frame(summarize(group_by(mcrpc_correlations, bin), 
  num_sig = sum(q_val < 0.05, na.rm = T), total = sum(!is.na(p_val))))

# Calculate proportion of significant correlations
mcrpc_all_correlations_binned$prop_sig = mcrpc_all_correlations_binned$num_sig/mcrpc_all_correlations_binned$total

# Create barplot of significant correlations per bin and the percentage of significant correlations covered by probes
mcrpc_prop_sig_bins_plot = ggplot(filter(mcrpc_all_correlations_binned, abs(bin) < 5000), aes(x = bin, y = prop_sig)) +
  geom_col(position = "identity", fill = "#2a5674") 
mcrpc_prop_sig_bins_plot = customize_ggplot_theme(mcrpc_prop_sig_bins_plot, xlab = "Distance to TSS", 
  ylab = "Proportion of CpG Sites\nDisplaying Significant Correlations", title = NULL)
mcrpc_prop_sig_bins_plot = mcrpc_prop_sig_bins_plot + scale_x_continuous(expand = c(0, 0), labels = scales::comma)
mcrpc_prop_sig_bins_plot
saveRDS(mcrpc_prop_sig_bins_plot, "mcrpc_prop_sig_bins_plot.rds")

### Create plots of standard deviations of CpG methylation

# Load CPGEA and MCRPC meth RSEs
cpgea_wgbs_hg38 = HDF5Array::loadHDF5SummarizedExperiment("../auxillary_data/methylation_data/cpgea_meth_rse")
mcrpc_wgbs_hg38 = HDF5Array::loadHDF5SummarizedExperiment("../auxillary_data/methylation_data/mcrpc_meth_rse")

# Load TSS and create TSS-proximal regions
tss_gr = readRDS("../auxillary_data/genomic_annotation/pc_tss_gr.rds")
tss_proximal_regions = methodical::expand_granges(tss_gr, 5000, 5000)

# Create 500 bp bins around TSS
bin_centers = seq(-4500, 4500, 500)
names(bin_centers) = as.character(bin_centers)
shiftStranded = function(x, shift=0L,...){
  GenomicRanges::shift(x ,shift=shift*ifelse('-'==strand(x),-1,1),...)}
tss_bins = lapply(bin_centers, function(x)
  shiftStranded(resize(tss_gr, 500, fix = "center"), x))

# Subset cpgea_wgbs_hg38 and mcrpc_wgbs_hg38 for CpGs overlapping TSS proximal regions
cpgea_wgbs_hg38_tss_proximal = subsetByOverlaps(cpgea_wgbs_hg38, reduce(tss_proximal_regions, ignore.strand = T))
mcrpc_wgbs_hg38_tss_proximal = subsetByOverlaps(mcrpc_wgbs_hg38, reduce(tss_proximal_regions, ignore.strand = T))

# Subset for normal samples
cpgea_wgbs_hg38_tss_proximal_normal = cpgea_wgbs_hg38_tss_proximal[, grep("N", colnames(cpgea_wgbs_hg38), value = T)]

# Create a GRanges with methylation standard deviations for all CpGs in TSS regions in normal samples. Took 6 minutes
system.time({cpg_sds_normal = DelayedMatrixStats::rowSds(assay(cpgea_wgbs_hg38_tss_proximal_normal), na.rm = T)})
names(cpg_sds_normal) = as.character(rowRanges(cpgea_wgbs_hg38_tss_proximal_normal))
cpg_sds_normal_gr = GRanges(names(cpg_sds_normal), values = unname(cpg_sds_normal))
saveRDS(cpg_sds_normal_gr, "cpg_sds_normal_gr.rds")

# Subset for tumour samples
cpgea_wgbs_hg38_tss_proximal_tumour = cpgea_wgbs_hg38_tss_proximal[, grep("T", colnames(cpgea_wgbs_hg38), value = T)]

# Create a GRanges with methylation standard deviations for all CpGs in TSS regions in tumour samples. . Took 6 minutes
system.time({cpg_sds_tumour = DelayedMatrixStats::rowSds(assay(cpgea_wgbs_hg38_tss_proximal_tumour), na.rm = T)})
names(cpg_sds_tumour) = as.character(rowRanges(cpgea_wgbs_hg38_tss_proximal_tumour))
cpg_sds_tumour_gr = GRanges(names(cpg_sds_tumour), values = unname(cpg_sds_tumour))
saveRDS(cpg_sds_tumour_gr, "cpg_sds_tumour_gr.rds")

# Create a GRanges with methylation standard deviations for all CpGs in TSS regions in metastasis samples
system.time({cpg_sds_metastasis = DelayedMatrixStats::rowSds(assay(mcrpc_wgbs_hg38_tss_proximal), na.rm = T)})
names(cpg_sds_metastasis) = as.character(rowRanges(mcrpc_wgbs_hg38_tss_proximal))
cpg_sds_metastasis_gr = GRanges(names(cpg_sds_metastasis), values = unname(cpg_sds_metastasis))
saveRDS(cpg_sds_metastasis_gr, "cpg_sds_metastasis_gr.rds")

# Load GRanges with standard deviations for CpG methylation
cpg_sds_normal_gr = readRDS("cpg_sds_normal_gr.rds")
cpg_sds_tumour_gr = readRDS("cpg_sds_tumour_gr.rds")
cpg_sds_metastasis_gr = readRDS("cpg_sds_metastasis_gr.rds")

# Create boxplots for SD of CpGs in bins around TSS for normal samples
cpg_sds_normal_bins = lapply(tss_bins, function(x)
  subsetByOverlaps(cpg_sds_normal_gr, x)$values)
cpg_sds_normal_bins = data.frame(
  bin = rep(names(cpg_sds_normal_bins), lengths(cpg_sds_normal_bins)),
  values = unlist(cpg_sds_normal_bins), row.names = NULL)
cpg_sds_normal_bins$bin = factor(cpg_sds_normal_bins$bin, levels = names(bin_centers))
cpg_sds_normal_boxplots = ggplot(cpg_sds_normal_bins, aes(x = as.factor(bin), y = values)) +
  geom_boxplot(outlier.shape = NA, fill = "#9E0142") 
cpg_sds_normal_boxplots = customize_ggplot_theme(cpg_sds_normal_boxplots, 
  ylab = "CpG Methylation Standard Deviation", xlab = "Distance to TSS (bp)",
   x_labels = c(rep("", 4), "-2,500", rep("", 4), "0", rep("", 4), "2,500", rep("", 4))) +
  scale_y_continuous(limits = c(0, 0.3), expand = c(0, 0))
cpg_sds_normal_boxplots
saveRDS(cpg_sds_normal_boxplots, "cpg_sds_normal_boxplots.rds")

# Create boxplots for SD of CpGs in bins around TSS for tumour samples
cpg_sds_tumour_bins = lapply(tss_bins, function(x)
  subsetByOverlaps(cpg_sds_tumour_gr, x)$values)
cpg_sds_tumour_bins = data.frame(
  bin = rep(names(cpg_sds_tumour_bins), lengths(cpg_sds_tumour_bins)),
  values = unlist(cpg_sds_tumour_bins), row.names = NULL)
cpg_sds_tumour_bins$bin = as.integer(cpg_sds_tumour_bins$bin)
cpg_sds_tumour_boxplots = ggplot(cpg_sds_tumour_bins, aes(x = as.factor(bin), y = values)) +
  geom_boxplot(outlier.shape = NA, fill = "#9E0142") 
cpg_sds_tumour_boxplots = customize_ggplot_theme(cpg_sds_tumour_boxplots, 
  ylab = "CpG Methylation Standard Deviation", xlab = "Distance to TSS (bp)", 
  x_labels = c(rep("", 4), "-2,500", rep("", 4), "0", rep("", 4), "2,500", rep("", 4))) +
  scale_y_continuous(limits = c(0, 0.4), expand = c(0, 0))
cpg_sds_tumour_boxplots
saveRDS(cpg_sds_tumour_boxplots, "cpg_sds_tumour_boxplots.rds")

# Create boxplots for SD of CpGs in bins around TSS for metastasis samples
cpg_sds_metastasis_bins = lapply(tss_bins, function(x)
  subsetByOverlaps(cpg_sds_metastasis_gr, x)$values)
cpg_sds_metastasis_bins = data.frame(
  bin = rep(names(cpg_sds_metastasis_bins), lengths(cpg_sds_metastasis_bins)),
  values = unlist(cpg_sds_metastasis_bins), row.names = NULL)
cpg_sds_metastasis_bins$bin = as.integer(cpg_sds_metastasis_bins$bin)
cpg_sds_metastasis_boxplots = ggplot(cpg_sds_metastasis_bins, aes(x = as.factor(bin), y = values)) +
  geom_boxplot(outlier.shape = NA, fill = "#9E0142") 
cpg_sds_metastasis_boxplots = customize_ggplot_theme(cpg_sds_metastasis_boxplots, 
  ylab = "CpG Methylation Standard Deviation", xlab = "Distance to TSS (bp)",
   x_labels = c(rep("", 4), "-2,500", rep("", 4), "0", rep("", 4), "2,500", rep("", 4))) +
  scale_y_continuous(limits = c(0, 0.55), expand = c(0, 0))
cpg_sds_metastasis_boxplots
saveRDS(cpg_sds_metastasis_boxplots, "cpg_sds_metastasis_boxplots.rds")

### Combine plots tom make figures

# Load plots for making combined plots
cpgs_per_bin_plot = readRDS("cpgs_per_bin_plot.rds")
combined_probes_per_bin_plot = readRDS("combined_probes_per_bin_plot.rds")
cpgea_normal_prop_sig_bins_plot = readRDS("cpgea_normal_prop_sig_bins_plot.rds")
cpgea_tumour_prop_sig_bins_plot = readRDS("cpgea_tumour_prop_sig_bins_plot.rds")
mcrpc_prop_sig_bins_plot = readRDS("mcrpc_prop_sig_bins_plot.rds")
cpg_sds_tumour_boxplots = readRDS("cpg_sds_tumour_boxplots.rds")
cpg_sds_normal_boxplots = readRDS("cpg_sds_normal_boxplots.rds")
cpg_sds_metastasis_boxplots = readRDS("cpg_sds_metastasis_boxplots.rds")

# Make combined plot with CpGs, probes and tumour plots
figure7_plot = cowplot::plot_grid(plotlist = 
    list(cpgs_per_bin_plot, combined_probes_per_bin_plot, cpgea_tumour_prop_sig_bins_plot, cpg_sds_tumour_boxplots), 
  nrow = 2, ncol = 2, align = "hv", labels = c("A", "B", "C", "D"), byrow = T)
ggsave(plot = figure7_plot, "../figures/figure7.pdf", width = 32, height = 18)

# Combine normal and metastases plots 
normal_and_metastases_plots = cowplot::plot_grid(plotlist = 
    list(cpgea_normal_prop_sig_bins_plot, cpg_sds_normal_boxplots,
      mcrpc_prop_sig_bins_plot, cpg_sds_metastasis_boxplots),
  nrow = 2, ncol = 2, align = "hv", labels = c("A", "B", "C", "D"), byrow = F)
normal_and_metastases_plots
ggsave(plot = normal_and_metastases_plots, "../figures/supp_normal_metastasis_probe_plot.pdf", width = 32, height = 18)

### Find proportion of TMRs overlapping Illumina probes

# Load TMRs for different datasets
cpgea_normal_tmrs = readRDS("../finding_tmrs/cpgea_normal_tmrs.rds")
cpgea_tumour_tmrs = readRDS("../finding_tmrs/cpgea_tumour_tmrs.rds")
mcrpc_tmrs = readRDS("../finding_tmrs/mcrpc_tmrs.rds")

# 30%, 40% and 46% of TMRs overlap Infinium 450K probes
length(subsetByOverlaps(cpgea_normal_tmrs, illumina_450k_probes_hg38))/length(cpgea_normal_tmrs)
length(subsetByOverlaps(cpgea_tumour_tmrs, illumina_450k_probes_hg38))/length(cpgea_tumour_tmrs)
length(subsetByOverlaps(mcrpc_tmrs, illumina_450k_probes_hg38))/length(mcrpc_tmrs)