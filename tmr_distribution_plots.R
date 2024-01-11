# Show spatial distribution of TMRs before and after removing repeats

# Load required packages
library(methodical)
library(dplyr)
library(plotR)
library(ggpubr)

# Get paths to all confident (those with at least 5 CpGs) TMR GRanges and put them in order
confident_tmr_granges_files = list.files("tmr_granges", full.names = T, pattern = "tmrs_5")
names(confident_tmr_granges_files) = basename(tools::file_path_sans_ext(confident_tmr_granges_files))
confident_tmr_granges_files = confident_tmr_granges_files[gtools::mixedorder(confident_tmr_granges_files)]

# Make a function which will bin relative TMRs
bin_relative_tmrs = function(tmr_file, width, transcripts_subset = NULL){
  
  # Load TMRs from file
  tmrs = readRDS(tmr_file)
  
  # If transcripts_subset provided, subet for TMRs associated with these transcripts
  if(!is.null(transcripts_subset)){
      tmrs = tmrs[tmrs$transcript_id %in% transcripts_subset]
  }
  
  # Convert the clusters to relative ranges
  relative_ranges_tmrs = 
    methodical::rangesRelativeToTSS(tmrs, tss_gr = GRanges(tmrs$tss_location))

  # Bin the relative ranges into 500 bp windows
  binned_ranges_tmrs = genomeTools::bin_relative_ranges(relative_ranges = relative_ranges_tmrs, bin_start = -width, 
    bin_end = width, bin_step = 500, category = tmrs$direction)
  
  # Convert to long format with separate rows for negative and positive TMRs
  binned_ranges_tmrs = reshape2::melt(binned_ranges_tmrs, id.vars = "bin_center")
  
  return(binned_ranges_tmrs)
  
}

# Bin relative TMRs for CPGEA Normal, CPGEA Tumour and MCRPC in 5KB
cpgea_normal_5kb_tmr_distributions_confident = bin_relative_tmrs(confident_tmr_granges_files["cpgea_normal_tmrs_5kb"], width = 4750)
cpgea_tumour_5kb_tmr_distributions_confident = bin_relative_tmrs(confident_tmr_granges_files["cpgea_tumour_tmrs_5kb"], width = 4750)
mcrpc_tmr_5kb_distributions_confident = bin_relative_tmrs(confident_tmr_granges_files["mcrpc_tmrs_5kb"], width = 4750)

# Combine TMR distributions_confident
combined_5kb_tmr_distributions_confident = bind_rows(
  cpgea_normal = cpgea_normal_5kb_tmr_distributions_confident, 
  cpgea_tumour = cpgea_tumour_5kb_tmr_distributions_confident,
  mcrpc = mcrpc_tmr_5kb_distributions_confident, .id = "dataset"
  )

# Make a plot of the TMR distrutions in the 3 data sets and save
tmrs_5kb_bins_plot_confident = ggplot(combined_5kb_tmr_distributions_confident, aes(x = bin_center, y = value, fill = variable)) +
  geom_col(position = "dodge")
tmrs_5kb_bins_plot_confident = customize_ggplot_theme(plot = tmrs_5kb_bins_plot_confident, base_theme = theme_bw(), xlab = "Distance to TSS (bp)", ylab = "Number of TMRs",
  title = NULL, fill_title = "TMR Direction", fill_colors = colour_list$purple_and_gold_light, axis_text_size = 14,
  fill_labels = c("Negative", "Positive"), scale_x = scale_x_continuous(breaks = seq(-4000, 4000, 2000), expand = c(0, 0), labels = scales::comma),
  facet = "dataset", facet_nrow = 1, facet_scales = "fixed", facet_labels = c("Normal Prostate", "Prostate Tumours", "Prostate Metastases"), strip_text_size = 18) + 
  theme(strip.background = element_blank(), plot.margin = margin(b = 0.5, unit = "cm")) +
  geom_vline(xintercept = 0, linetype = "dotted")
tmrs_5kb_bins_plot_confident + labs(fill = "Cluster Direction", y = "Number of Clusters")
ggsave(plot = tmrs_5kb_bins_plot_confident , "tmr_distribution_plots/confident_tmrs_5kb_bins_plot.pdf", width = 32, height  = 9)
ggsave(plot = tmrs_5kb_bins_plot_confident + labs(fill = "Cluster Direction", y = "Number of Clusters") , 
  "../presentation_figures/confident_tmrs_5kb_bins_plot.pdf", width = 32, height  = 9)

# Bin relative TMRs for CPGEA Normal, CPGEA Tumour and MCRPC in 50KB
cpgea_normal_50kb_tmr_distributions_confident = bin_relative_tmrs(confident_tmr_granges_files["cpgea_normal_tmrs_50kb"], width = 50000-250)
cpgea_tumour_50kb_tmr_distributions_confident = bin_relative_tmrs(confident_tmr_granges_files["cpgea_tumour_tmrs_50kb"], width = 50000-250)
mcrpc_tmr_50kb_distributions_confident = bin_relative_tmrs(confident_tmr_granges_files["mcrpc_tmrs_50kb"], width = 50000-250)

# Combine TMR distributions
combined_50kb_tmr_distributions_confident = bind_rows(
  cpgea_normal = cpgea_normal_50kb_tmr_distributions_confident, 
  cpgea_tumour = cpgea_tumour_50kb_tmr_distributions_confident,
  mcrpc = mcrpc_tmr_50kb_distributions_confident, .id = "dataset"
  )

# Make a plot of the TMR distrutions in the 3 data sets and save
tmrs_50kb_bins_plot_confident = ggplot(combined_50kb_tmr_distributions_confident, aes(x = bin_center, y = value, fill = variable)) +
  geom_col(position = "dodge", color = NA)
tmrs_50kb_bins_plot_confident = customize_ggplot_theme(plot = tmrs_50kb_bins_plot_confident, base_theme = theme_bw(), xlab = "Distance to TSS (bp)", ylab = "Number of TMRs",
  title = NULL, fill_title = "TMR Direction", fill_colors = colour_list$purple_and_gold_light, axis_text_size = 14,
  fill_labels = c("Negative", "Positive"), scale_x = scale_x_continuous(breaks = seq(-40000, 40000, 20000), expand = c(0, 0), labels = scales::comma),
  facet = "dataset", facet_nrow = 1, facet_scales = "fixed", facet_labels = c("Normal Prostate", "Prostate Tumours", "Prostate Metastases"), strip_text_size = 18) + 
  theme(strip.background = element_blank(), plot.margin = margin(b = 0.5, unit = "cm")) +
  geom_vline(xintercept = 0, linetype = "dotted")
ggsave(plot = tmrs_50kb_bins_plot_confident + labs(fill = "Cluster Direction", y = "Number of Clusters"), 
  "../presentation_figures/confident_tmrs_50kb_bins_plot.pdf", width = 32, height  = 9)
ggsave(plot = tmrs_50kb_bins_plot_confident , "tmr_distribution_plots/confident_tmrs_50kb_bins_plot.pdf", width = 32, height  = 9)

# Bin relative TMRs for CPGEA Normal, CPGEA Tumour and MCRPC in 500KB
cpgea_normal_500kb_tmr_distributions_confident = bin_relative_tmrs(confident_tmr_granges_files["cpgea_normal_tmrs_500kb"], width = 500000-250)
cpgea_tumour_500kb_tmr_distributions_confident = bin_relative_tmrs(confident_tmr_granges_files["cpgea_tumour_tmrs_500kb"], width = 500000-250)
mcrpc_tmr_500kb_distributions_confident = bin_relative_tmrs(confident_tmr_granges_files["mcrpc_tmrs_500kb"], width = 500000-250)

# Combine TMR distributions
combined_500kb_tmr_distributions_confident = bind_rows(
  cpgea_normal = cpgea_normal_500kb_tmr_distributions_confident, 
  cpgea_tumour = cpgea_tumour_500kb_tmr_distributions_confident,
  mcrpc = mcrpc_tmr_500kb_distributions_confident, .id = "dataset"
  )

# Make a plot of the TMR distrutions in the 3 data sets and save
tmrs_500kb_bins_plot_confident = ggplot(combined_500kb_tmr_distributions_confident, aes(x = bin_center, y = value, fill = variable)) +
  geom_col(position = "dodge", color = NA)
tmrs_500kb_bins_plot_confident = customize_ggplot_theme(plot = tmrs_500kb_bins_plot_confident, base_theme = theme_bw(), xlab = "Distance to TSS (bp)", ylab = "Number of TMRs",
  title = NULL, fill_title = "TMR Direction", fill_colors = colour_list$purple_and_gold_light, axis_text_size = 14,
  fill_labels = c("Negative", "Positive"), scale_x = scale_x_continuous(breaks = seq(-400000, 400000, 200000), expand = c(0, 0), labels = scales::comma),
  facet = "dataset", facet_nrow = 1, facet_scales = "fixed", facet_labels = c("Normal Prostate", "Prostate Tumours", "Prostate Metastases"), strip_text_size = 18) + 
  theme(strip.background = element_blank(), plot.margin = margin(b = 0.5, unit = "cm")) +
  geom_vline(xintercept = 0, linetype = "dotted", linewidth = 0.5, alpha = 0.25)
ggsave(plot = tmrs_500kb_bins_plot_confident , "tmr_distribution_plots/confident_tmrs_500kb_bins_plot.pdf", width = 32, height  = 9)

# Get paths to all final (those without repeats) TMR GRanges and put them in order
final_tmr_granges_files = list.files("final_tmr_granges", full.names = T)
names(final_tmr_granges_files) = basename(tools::file_path_sans_ext(final_tmr_granges_files))
final_tmr_granges_files = final_tmr_granges_files[gtools::mixedorder(final_tmr_granges_files)]

# Bin relative TMRs for CPGEA Normal, CPGEA Tumour and MCRPC in 5KB
cpgea_normal_5kb_tmr_distributions_final = bin_relative_tmrs(final_tmr_granges_files["cpgea_normal_tmrs_5kb"], width = 4750)
cpgea_tumour_5kb_tmr_distributions_final = bin_relative_tmrs(final_tmr_granges_files["cpgea_tumour_tmrs_5kb"], width = 4750)
mcrpc_tmr_5kb_distributions_final = bin_relative_tmrs(final_tmr_granges_files["mcrpc_tmrs_5kb"], width = 4750)

# Combine TMR distributions_final
combined_5kb_tmr_distributions_final = bind_rows(
  cpgea_normal = cpgea_normal_5kb_tmr_distributions_final, 
  cpgea_tumour = cpgea_tumour_5kb_tmr_distributions_final,
  mcrpc = mcrpc_tmr_5kb_distributions_final, .id = "dataset"
  )

# Make a plot of the TMR distrutions in the 3 data sets and save
tmrs_5kb_bins_plot_final = ggplot(combined_5kb_tmr_distributions_final, aes(x = bin_center, y = value, fill = variable)) +
  geom_col(position = "dodge")
tmrs_5kb_bins_plot_final = customize_ggplot_theme(plot = tmrs_5kb_bins_plot_final, base_theme = theme_bw(), xlab = "Distance to TSS (bp)", ylab = "Number of TMRs",
  title = NULL, fill_title = "TMR Direction", fill_colors = colour_list$purple_and_gold_light, axis_text_size = 14,
  fill_labels = c("Negative", "Positive"), scale_x = scale_x_continuous(breaks = seq(-4000, 4000, 2000), expand = c(0, 0), labels = scales::comma),
  facet = "dataset", facet_nrow = 1, facet_scales = "fixed", facet_labels = c("Normal Prostate", "Prostate Tumours", "Prostate Metastases"), strip_text_size = 18) + 
  theme(strip.background = element_blank(), plot.margin = margin(b = 0.5, unit = "cm")) +
  geom_vline(xintercept = 0, linetype = "dotted")
ggsave(plot = tmrs_5kb_bins_plot_final + labs(fill = "Cluster Direction", y = "Number of Clusters"), 
  "../presentation_figures//final_tmrs_5kb_bins_plot.pdf", width = 24, height  = 6.75)
ggsave(plot = tmrs_5kb_bins_plot_final , "tmr_distribution_plots/final_tmrs_5kb_bins_plot.pdf", width = 32, height  = 9)
saveRDS(tmrs_5kb_bins_plot_final, "tmr_distribution_plots/tmrs_5kb_bins_plot_final.rds")

# Bin relative TMRs for CPGEA Normal, CPGEA Tumour and MCRPC in 50KB
cpgea_normal_50kb_tmr_distributions_final = bin_relative_tmrs(final_tmr_granges_files["cpgea_normal_tmrs_50kb"], width = 50000-250)
cpgea_tumour_50kb_tmr_distributions_final = bin_relative_tmrs(final_tmr_granges_files["cpgea_tumour_tmrs_50kb"], width = 50000-250)
mcrpc_tmr_50kb_distributions_final = bin_relative_tmrs(final_tmr_granges_files["mcrpc_tmrs_50kb"], width = 50000-250)

# Combine TMR distributions
combined_50kb_tmr_distributions_final = bind_rows(
  cpgea_normal = cpgea_normal_50kb_tmr_distributions_final, 
  cpgea_tumour = cpgea_tumour_50kb_tmr_distributions_final,
  mcrpc = mcrpc_tmr_50kb_distributions_final, .id = "dataset"
  )

# Make a plot of the TMR distrutions in the 3 data sets and save
tmrs_50kb_bins_plot_final = ggplot(combined_50kb_tmr_distributions_final, aes(x = bin_center, y = value, fill = variable)) +
  geom_col(position = "dodge", color = NA)
tmrs_50kb_bins_plot_final = customize_ggplot_theme(plot = tmrs_50kb_bins_plot_final, base_theme = theme_bw(), xlab = "Distance to TSS (bp)", ylab = "Number of TMRs",
  title = NULL, fill_title = "TMR Direction", fill_colors = colour_list$purple_and_gold_light, axis_text_size = 14,
  fill_labels = c("Negative", "Positive"), scale_x = scale_x_continuous(breaks = seq(-40000, 40000, 20000), expand = c(0, 0), labels = scales::comma),
  facet = "dataset", facet_nrow = 1, facet_scales = "fixed", facet_labels = c("Normal Prostate", "Prostate Tumours", "Prostate Metastases"), strip_text_size = 18) + 
  theme(strip.background = element_blank(), plot.margin = margin(b = 0.5, unit = "cm")) +
  geom_vline(xintercept = 0, linetype = "dotted") +
  geom_vline(xintercept = -5000, linetype = "dotted") +
  geom_vline(xintercept = 5000, linetype = "dotted")
ggsave(plot = tmrs_50kb_bins_plot_final + labs(fill = "Cluster Direction", y = "Number of Clusters"), 
  "../presentation_figures//tmrs_50kb_bins_plot_final.pdf", width = 24, height  = 6.75)
ggsave(plot = tmrs_50kb_bins_plot_final , "tmr_distribution_plots/final_tmrs_50kb_bins_plot.pdf", width = 32, height  = 9)

# Bin relative TMRs for CPGEA Normal, CPGEA Tumour and MCRPC in 500KB
cpgea_normal_500kb_tmr_distributions_final = bin_relative_tmrs(final_tmr_granges_files["cpgea_normal_tmrs_500kb"], width = 500000-250)
cpgea_tumour_500kb_tmr_distributions_final = bin_relative_tmrs(final_tmr_granges_files["cpgea_tumour_tmrs_500kb"], width = 500000-250)
mcrpc_tmr_500kb_distributions_final = bin_relative_tmrs(final_tmr_granges_files["mcrpc_tmrs_500kb"], width = 500000-250)

# Combine TMR distributions
combined_500kb_tmr_distributions_final = bind_rows(
  cpgea_normal = cpgea_normal_500kb_tmr_distributions_final, 
  cpgea_tumour = cpgea_tumour_500kb_tmr_distributions_final,
  mcrpc = mcrpc_tmr_500kb_distributions_final, .id = "dataset"
  )

# Make a plot of the TMR distrutions in the 3 data sets and save
tmrs_500kb_bins_plot_final = ggplot(combined_500kb_tmr_distributions_final, aes(x = bin_center, y = value, fill = variable)) +
  geom_col(position = "dodge", color = NA)
tmrs_500kb_bins_plot_final = customize_ggplot_theme(plot = tmrs_500kb_bins_plot_final, base_theme = theme_bw(), xlab = "Distance to TSS (bp)", ylab = "Number of TMRs",
  title = NULL, fill_title = "TMR Direction", fill_colors = colour_list$purple_and_gold_light, axis_text_size = 14,
  fill_labels = c("Negative", "Positive"), scale_x = scale_x_continuous(breaks = seq(-400000, 400000, 200000), expand = c(0, 0), labels = scales::comma),
  facet = "dataset", facet_nrow = 1, facet_scales = "fixed", facet_labels = c("Normal Prostate", "Prostate Tumours", "Prostate Metastases"), strip_text_size = 18) + 
  theme(strip.background = element_blank(), plot.margin = margin(b = 0.5, unit = "cm")) +
  geom_vline(xintercept = 0, linetype = "dotted", linewidth = 0.5, alpha = 0.25)
ggsave(plot = tmrs_500kb_bins_plot_final , "tmr_distribution_plots/final_tmrs_500kb_bins_plot.pdf", width = 32, height  = 9)

# Add lines showing 5 KB on tmrs_50kb_bins_plot_final
tmrs_50kb_bins_plot_final_with_5kb_lines = tmrs_50kb_bins_plot_final + 
  geom_vline(xintercept = -5000, linetype = "dashed") + 
  geom_vline(xintercept = 5000, linetype = "dashed")
ggsave(plot = tmrs_50kb_bins_plot_final_with_5kb_lines , "tmr_distribution_plots/final_tmrs_50kb_bins_plot_final_with_5kb_lines.pdf", width = 32, height  = 9)

# Combine all TMR plots with repeats into a single plot
repeat_distribution_plot_list = list(tmrs_5kb_bins_plot_confident, tmrs_50kb_bins_plot_confident, tmrs_500kb_bins_plot_confident)
repeat_distribution_plots = ggpubr::ggarrange(plotlist = repeat_distribution_plot_list, nrow = 3, labels = c("A", "B", "C"))
ggsave(plot = repeat_distribution_plots , "tmr_distribution_plots/repeat_distribution_plots.pdf", width = 32, height  = 27)

# Combine all TMR plots without repeats into a single plot
no_repeat_distribution_plot_list = list(tmrs_5kb_bins_plot_final, tmrs_50kb_bins_plot_final, tmrs_500kb_bins_plot_final)
no_repeat_distribution_plots = ggarrange(plotlist = no_repeat_distribution_plot_list, nrow = 3, labels = c("A", "B", "C"))
ggsave(plot = no_repeat_distribution_plots , "tmr_distribution_plots/no_repeat_distribution_plots.pdf", width = 32, height  = 27)

# Create a theme for use with poster
poster_theme = theme(plot.title = element_text(hjust = 0.5, size = 30), 
	axis.title = element_text(size = 26), axis.text = element_text(size = 24), strip.text.x = element_text(size = 30),
	legend.text = element_text(size = 24), legend.title = element_text(size = 26), legend.key.size = unit(2, "cm")) 

# Make a plot with 5KB and 50KB plots for poster
no_repeat_distribution_plots_poster = ggarrange(plotlist = 
    list(tmrs_5kb_bins_plot_final + poster_theme, tmrs_50kb_bins_plot_final_with_5kb_lines + poster_theme), nrow = 2, labels = c("A", "B"))
ggsave(plot = no_repeat_distribution_plots_poster, "~/promoter_project/presentation_figures/no_repeat_distribution_plots_poster.pdf", width = 32, height  = 18)

### Canonical transcripts evaluation

# Get list of canonical transcripts from Ensembl
canonical_transcripts = readRDS("~/genomes/gencode/gencode_granges/ensembl_canonical_transcripts.rds")

# Bin relative TMRs for CPGEA Normal, CPGEA Tumour and MCRPC in 5KB
cpgea_normal_5kb_tmr_distributions_canonical = bin_relative_tmrs(tmr_file = final_tmr_granges_files["cpgea_normal_tmrs_5kb"], width = 4750, transcripts_subset = canonical_transcripts)
cpgea_tumour_5kb_tmr_distributions_canonical = bin_relative_tmrs(final_tmr_granges_files["cpgea_tumour_tmrs_5kb"], width = 4750, transcripts_subset = canonical_transcripts)
mcrpc_tmr_5kb_distributions_canonical = bin_relative_tmrs(final_tmr_granges_files["mcrpc_tmrs_5kb"], width = 4750, transcripts_subset = canonical_transcripts)

# Combine TMR distributions_canonical
combined_5kb_tmr_distributions_canonical = bind_rows(
  cpgea_normal = cpgea_normal_5kb_tmr_distributions_canonical, 
  cpgea_tumour = cpgea_tumour_5kb_tmr_distributions_canonical,
  mcrpc = mcrpc_tmr_5kb_distributions_canonical, .id = "dataset"
  )

# Make a plot of the TMR distrutions in the 3 data sets and save
tmrs_5kb_bins_plot_canonical = ggplot(combined_5kb_tmr_distributions_canonical, aes(x = bin_center, y = value, fill = variable)) +
  geom_col(position = "dodge")
tmrs_5kb_bins_plot_canonical = customize_ggplot_theme(plot = tmrs_5kb_bins_plot_canonical, base_theme = theme_bw(), xlab = "Distance to TSS (bp)", ylab = "Number of TMRs",
  title = NULL, fill_title = "TMR Direction", fill_colors = colour_list$purple_and_gold_light, axis_text_size = 14,
  fill_labels = c("Negative", "Positive"), scale_x = scale_x_continuous(breaks = seq(-4000, 4000, 2000), expand = c(0, 0), labels = scales::comma),
  facet = "dataset", facet_nrow = 1, facet_scales = "fixed", facet_labels = c("Normal Prostate", "Prostate Tumours", "Prostate Metastases"), strip_text_size = 18) + 
  theme(strip.background = element_blank(), plot.margin = margin(b = 0.5, unit = "cm")) +
  geom_vline(xintercept = 0, linetype = "dotted")
ggsave(plot = tmrs_5kb_bins_plot_canonical , "tmr_distribution_plots/canonical_tmrs_5kb_bins_plot.pdf", width = 32, height  = 9)

### Presentation plots
# Create a theme for use with poster
poster_theme = theme(plot.title = element_text(hjust = 0.5, size = 30), 
	axis.title = element_text(size = 26), axis.text = element_text(size = 22), strip.text = element_text(size = 30),
	legend.text = element_text(size = 24), legend.title = element_text(size = 26), legend.key.size = unit(2, "cm")) 

repeat_distribution_plots_poster = ggarrange(plotlist = lapply(repeat_distribution_plot_list, function(x) x + poster_theme), nrow = 3, labels = c("A", "B", "C"), 
  font.label = list(size = 26))
ggsave(plot = repeat_distribution_plots_poster, "~/promoter_project/martin_plots/repeat_distribution_plots.pdf", width = 32, height  = 27)

### Rhabdoid plots
rhabdoid_tmr_distribution = bin_relative_tmrs(final_tmr_granges_files["rhabdoid_tmrs_5kb"], width = 4750)
rhabdoid_tmrs_bin_plot = ggplot(rhabdoid_tmr_distribution, aes(x = bin_center, y = value, fill = variable)) +
  geom_col(position = "dodge")
customize_ggplot_theme(plot = rhabdoid_tmrs_bin_plot, base_theme = theme_bw(), xlab = "Distance to TSS (bp)", ylab = "Number of TMRs",
  title = NULL, fill_title = "TMR Direction", fill_colors = colour_list$purple_and_gold_light, axis_text_size = 14,
  fill_labels = c("Negative", "Positive"), scale_x = scale_x_continuous(breaks = seq(-4000, 4000, 2000), expand = c(0, 0), labels = scales::comma)) + 
  theme(strip.background = element_blank(), plot.margin = margin(b = 0.5, unit = "cm")) +
  geom_vline(xintercept = 0, linetype = "dotted", linewidth = 0.5, alpha = 0.25)