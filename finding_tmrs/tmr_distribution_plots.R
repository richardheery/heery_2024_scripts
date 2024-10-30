# Create TMR distribution plots

# Load required packages
library(GenomicRanges)
library(dplyr)
library(ggpubr)
library(doParallel)
source("../auxillary_scripts/plotting_functions.R")
source("../auxillary_scripts/tmr_plot_functions.R")

# Load TMRs filtered for CAGE supported TSS not overlapping repeats
cpgea_normal_tmrs = readRDS("cpgea_normal_tmrs.rds")
cpgea_tumour_tmrs = readRDS("cpgea_tumour_tmrs.rds")
mcrpc_tmrs = readRDS("mcrpc_tmrs.rds")

# Bin relative TMRs for CPGEA Normal, CPGEA Tumour and MCRPC in 5KB
cpgea_normal_5kb_tmr_distributions = bin_relative_tmrs(cpgea_normal_tmrs, width = 4750)
cpgea_tumour_5kb_tmr_distributions = bin_relative_tmrs(cpgea_tumour_tmrs, width = 4750)
mcrpc_tmr_5kb_distributions = bin_relative_tmrs(mcrpc_tmrs, width = 4750)

# Combine TMR distributions
combined_5kb_tmr_distributions = bind_rows(
  cpgea_normal = cpgea_normal_5kb_tmr_distributions, 
  cpgea_tumour = cpgea_tumour_5kb_tmr_distributions,
  mcrpc = mcrpc_tmr_5kb_distributions, .id = "dataset"
  )

# Make a plot of the TMR distrutions in the 3 data sets and save
tmrs_5kb_bins_plot = ggplot(combined_5kb_tmr_distributions, aes(x = bin_center, y = value, fill = variable)) +
  geom_col(position = "dodge")
tmrs_5kb_bins_plot = customize_ggplot_theme(plot = tmrs_5kb_bins_plot, base_theme = theme_bw(), xlab = "Distance to TSS (bp)", ylab = "Number of TMRs",
  title = NULL, fill_title = "TMR Direction", fill_colors = colour_list$purple_and_gold_light, axis_text_size = 14,
  fill_labels = c("Negative", "Positive"), scale_x = scale_x_continuous(breaks = seq(-4000, 4000, 2000), expand = c(0, 0), labels = scales::comma),
  facet = "dataset", facet_nrow = 1, facet_scales = "fixed", facet_labels = c("Normal Prostate", "Prostate Tumours", "Prostate Metastases"), strip_text_size = 18) + 
  theme(strip.background = element_blank(), plot.margin = margin(b = 0.5, unit = "cm")) +
  geom_vline(xintercept = 0, linetype = "dotted")
ggsave(plot = tmrs_5kb_bins_plot , "../figures/supplementary_figure8A.pdf",  width = 27, height = 9)

### Make introns and exons plots

tmr_list = list(
  "Normal Prostate" = cpgea_normal_tmrs,
  "Prostate Tumours" = cpgea_tumour_tmrs,
  "Prostate Metastases" = mcrpc_tmrs)

# Load introns and exons for transcripts as a GRanges list
tss_grl = readRDS("../auxillary_data/genomic_annotation/pc_transcripts_exons_and_introns_grl.rds")
tss_grl_expanded = expand_transcripts(tss_grl, seq(500, 2000, 500), seq(500, 2000, 500))

# Create a vector with the regions to plot
exons = paste0("exon_", 1:10)
introns = paste0("intron_", 1:10)
exons_introns = c(rbind(exons, introns))
exons_introns = exons_introns[-length(exons_introns)]
regions = c(paste0("TSS-", rev(seq(500, 2000, 500))), exons_introns, paste0("TES+", seq(500, 2000, 500)))

# Create absolute count and normalized counts for CPGEA normal samples
cpgea_normal_tmrs_introns_exons_plot = plot_tmr_regions(tmrs = cpgea_normal_tmrs, 
  transcript_regions_gr = tss_grl_expanded, regions_filter = regions, title = "Distribution of TMRs in Normal Prostate", normalize = F)
cpgea_normal_tmrs_introns_exons_plot_normalized = plot_tmr_regions(tmrs = cpgea_normal_tmrs, 
  transcript_regions_gr = tss_grl_expanded, regions_filter = regions, title = "Distribution of TMRs in Normal Prostate", normalize = T)
ggsave(plot = cpgea_normal_tmrs_introns_exons_plot, "../figures/supplementary_figure9.pdf", width = 16, height = 9)
ggsave(plot = ggarrange(cpgea_normal_tmrs_introns_exons_plot_normalized, labels = "D"), 
  "../figures/figure4D.pdf", width = 20.57, height = 11.57)

# Create absolute count and normalized counts for CPGEA tumour samples
cpgea_tumour_tmrs_introns_exons_plot = plot_tmr_regions(tmrs = cpgea_tumour_tmrs, 
  transcript_regions_gr = tss_grl_expanded, regions_filter = regions, title = NULL, normalize = F)
cpgea_tumour_tmrs_introns_exons_plot_normalized = plot_tmr_regions(tmrs = cpgea_tumour_tmrs, 
  transcript_regions_gr = tss_grl_expanded, regions_filter = regions, title = NULL, normalize = T)

# Create absolute count and normalized counts for MCRPC samples
mcrpc_tmrs_introns_exons_plot = plot_tmr_regions(tmrs = mcrpc_tmrs, 
  transcript_regions_gr = tss_grl_expanded, regions_filter = regions, title = NULL, normalize = F)
mcrpc_tmrs_introns_exons_plot_normalized = plot_tmr_regions(tmrs = mcrpc_tmrs, 
  transcript_regions_gr = tss_grl_expanded, regions_filter = regions, title = NULL, normalize = T)

# Combine tumour and metastases plots
tumour_and_metastases_plots = list(cpgea_tumour_tmrs_introns_exons_plot, cpgea_tumour_tmrs_introns_exons_plot_normalized,
  mcrpc_tmrs_introns_exons_plot, mcrpc_tmrs_introns_exons_plot_normalized)
combined_tumour_metastases_plots = ggarrange(plotlist = tumour_and_metastases_plots, nrow = 2, ncol = 2, labels = LETTERS[1:4], common.legend = T, legend = "right")
combined_tumour_metastases_plots
ggsave(plot = combined_tumour_metastases_plots, "../figures/supplementary_figure10.pdf", width = 32, height = 18)

# Make distribution plots just for MANE transcripts
mane_transcripts = readRDS("../auxillary_data/genomic_annotation/mane_transcript_ids.rds")
cpgea_normal_mane_tmrs_introns_exons_plot_normalized = plot_tmr_regions(tmrs = cpgea_normal_tmrs[cpgea_normal_tmrs$ID %in% mane_transcripts], 
  transcript_regions_gr = tss_grl_expanded, regions_filter = regions, title = "Distribution of TMRs for MANE Transcripts in Normal Prostate", normalize = T)
ggsave(plot = cpgea_normal_mane_tmrs_introns_exons_plot_normalized, "../figures/supplementary_figure12.pdf", width = 16, height = 9)
