# Create plots of TMR stats

# Load required packages
library(plotR)
library(genomeTools)

# Load TMRs
cpgea_normal_tmrs = readRDS("final_tmr_granges/cpgea_normal_tmrs_5kb.rds")
cpgea_tumour_tmrs = readRDS("final_tmr_granges/cpgea_tumour_tmrs_5kb.rds")
mcrpc_tmrs = readRDS("final_tmr_granges/mcrpc_tmrs_5kb.rds")

# Create a list with the different TMR types
tmr_list = list(
  cpgea_normal = cpgea_normal_tmrs,
  cpgea_tumour = cpgea_tumour_tmrs,
  mcrpc = mcrpc_tmrs
)

# Get transcripts and genes associated with TMRs from each dataset
tmr_transcripts = lapply(tmr_list, function(x) unique(x$transcript_id))
tmr_genes = lapply(tmr_list, function(x) unique(x$gene_name))
lengths(tmr_transcripts)
lengths(tmr_genes)

# Get list of TMRs
tmr_list = readRDS("final_tmr_granges/final_tmr_5kb_list.rds")

# Create a data.frame summarizing the number of TMRs, transcripts and genes associated with each TMR group
tmr_stats = data.frame(
  dataset = gsub("_negative|_positive", "", names(tmr_list)),
  direction = stringr::str_to_title(gsub(".*_", "", names(tmr_list))),
  tmr_count = lengths(tmr_list),
  transcript_count = sapply(tmr_list, function(x) length(unique(x$transcript_id))),
  gene_count = sapply(tmr_list, function(x) length(unique(x$gene_name))),
  row.names = NULL
)

# Put dataset levels in correct order
tmr_stats$dataset = factor(tmr_stats$dataset, unique(tmr_stats$dataset))

# Convert tmr_stats into long format
tmr_stats = tidyr::pivot_longer(tmr_stats, cols = c("tmr_count", "transcript_count", "gene_count"))
tmr_stats$name = factor(tmr_stats$name, levels = c("tmr_count", "transcript_count", "gene_count"))

# Create barplots with the number of TMRs, TMR-associated transcripts and TMR-associated genes for each dataset
tmr_stats_barplot = ggplot(tmr_stats, aes(x = dataset, y = value, fill = direction)) +
  geom_col(position = "dodge", color = "black")
tmr_stats_barplot = customize_ggplot_theme(tmr_stats_barplot, 
  title = NULL, xlab = "Dataset", y = "Count", fill_title = "TMR Direction",
  fill_colors = colour_list$purple_and_gold_light, x_labels = c("Normal Prostate", "Prostate Tumours", "Prostate Metastases"), 
  facet = "name", facet_labels = c("Number of TMRs", "Number of TMR-Associated Transcripts", "Number of TMR-Associated Genes"), 
  facet_scales = "fixed", x_labels_angle = 30) + theme(strip.background = element_blank())
ggsave(plot = tmr_stats_barplot, filename = "tmr_stats_plots/tmr_stats_barplot.pdf", width = 27, height = 9)

# Calculate the proportion overlap between TMR groups
tmr_overlaps = calculate_regions_overlap_list(tmr_list, ignore.strand = T, overlap_threshold = 0.25)

# Get overlap proportions for transcripts and genes from different TMR groups
transcript_overlaps = intersectR::intersect_lengths_all_pairwise(lapply(tmr_list, function(x) x$transcript_id), proportion = T)
gene_overlaps = intersectR::intersect_lengths_all_pairwise(lapply(tmr_list, function(x) x$gene_name), proportion = T)

# Create row and column labels for heatmaps 
labels = c("Normal Prostate TMRs -", "Normal Prostate TMRs +", "Prostate Tumour TMRs -", 
  "Prostate Tumour TMRs +", "Prostate Metastasis TMRs -", "Prostate Metastasis TMRs +")

# Create heatmaps of TMR overlaps, for Jaccard indices for transcripts and genes
tmr_overlaps_plot = plotR::heatmap_without_clustering(tmr_overlaps, row_labels = labels, col_labels = labels, 
  filename = NA, title = "Relative Overlap of Genomic Regions\nCovered by TMRs", title_size = 15, return_ggplot = T)
transcript_overlaps_plot = plotR::heatmap_without_clustering(transcript_overlaps, row_labels = labels, col_labels = labels, 
  filename = NA, title = "Relative Overlap of Transcripts\nAssociated with TMRs", title_size = 15, return_ggplot = T)
gene_overlaps_plot = plotR::heatmap_without_clustering(gene_overlaps, row_labels = labels, col_labels = labels, 
  filename = NA, title = "Relative Overlap of Genes\nAssociated with TMRs", title_size = 15, return_ggplot = T)

# Combine heatmaps and save
tmr_heatmap_list = list(tmr_overlaps_plot, transcript_overlaps_plot, gene_overlaps_plot)
combined_tmr_heatmaps = ggpubr::ggarrange(plotlist = tmr_heatmap_list, nrow = 1, ncol = 3)
ggsave(plot = combined_tmr_heatmaps, filename = "tmr_stats_plots/combined_tmr_heatmaps.pdf", width = 27, height = 9, bg = "white")
ggsave(plot = tmr_overlaps_plot, filename = "~/promoter_project/presentation_figures/tmr_overlaps_plot.pdf", width = 9, height = 9)

# Create plots for presentation and save
tmr_count_barplot = ggplot(filter(tmr_stats, name == "tmr_count"), aes(x = dataset, y = value, fill = direction)) +
  geom_col(position = "dodge", color = "black")
tmr_count_barplot = customize_ggplot_theme(tmr_count_barplot , 
  title = NULL, xlab = "Dataset", y = "Number of TMRs", fill_title = "TMR Direction",
  fill_colors = colour_list$purple_and_gold_light, x_labels = c("Normal Prostate", "Prostate Tumours", "Prostate Metastases"))
ggsave(plot = tmr_count_barplot, filename = "~/promoter_project/presentation_figures/tmr_count_barplot.pdf", width = 16, height = 9)
ggsave(plot = tmr_overlaps_plot, filename = "~/promoter_project/presentation_figures/tmr_overlaps_plot.pdf", width = 9, height = 9)
