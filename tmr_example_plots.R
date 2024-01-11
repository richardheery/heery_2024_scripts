# Create plot of Methodical workflow and show examples of some TMRs

# Load required packages
library(methodical)
library(ggpubr)

# Load CPGEA normal correlation results
cpgea_normal_correlation_results = readRDS("~/promoter_project/meth_transcript_correlations/meth_transcript_cor_lists/transcript_meth_cors_cpgea_normal_samples_5kb.rds")

# Get GRanges for protein-coding TSS sites
tss_gr = readRDS("/media/rich/Elements/genomes/gencode/gencode_granges/pcg_transcript_tss_ranges_gencode_v38.rds")

# Gt GRanges for TMRs for normal prostate samples
cpgea_normal_tmrs = readRDS("~/promoter_project/final_tmrs/final_tmr_granges/cpgea_normal_tmrs_5kb.rds")

### Show TMR workflow using TUBB6 (ENST00000591909) as an example

# Get CpG correlation values for TUBB6
tubb6_cpg_correlations = cpgea_normal_correlation_results[["ENST00000591909"]]
tubb6_cpg_correlations  = dplyr::rename(tubb6_cpg_correlations, meth_site = "cpg_name")
#tubb6_cpg_correlations = tibble::column_to_rownames(tubb6_cpg_correlations, "meth_site")

# Get TUBB6 TSS site
tubb6_tss = plyranges::filter(tss_gr, transcript_id == "ENST00000591909")

# Get TMRs for TUBB6
tubb6_tmrs = plyranges::filter(cpgea_normal_tmrs, transcript_id == "ENST00000591909")

# Create TUBB6 CpG correlation plot
tubb6_cpg_correlation_plot = plotMethSiteCorCoefs(meth_site_cor_values = tubb6_cpg_correlations, reference_tss = tubb6_tss, 
  value_colours = "set2", xlabel = "Distance to TSS (bp)", ylabel = "DNA Methylation-Transcription Correlation") +
  scale_x_continuous(limits = c(-5000, 5000), breaks = seq(-4000, 4000, 2000), exp= c(0.005, 0.005), labels = scales::comma) +
  geom_hline(yintercept = 0, linetype = "dashed") + 
  labs(title = "1: Calculate CpG methylation-transcription correlation values")
ggsave(plot = tubb6_cpg_correlation_plot, "~/promoter_project/presentation_figures/methodical_part1.pdf", width = 16, height = 9)

# Extract the colours for the points
plot_colours = ggplot_build(tubb6_cpg_correlation_plot)$data[[2]][["fill"]]

# Create TUBB6 methodical scores plot
tubb6_methodical_scores_plot = plotMethodicalScores(meth_site_values = tubb6_cpg_correlations, reference_tss = tubb6_tss, xlabel = "Distance to TSS (bp)", 
  p_value_threshold = NULL, smooth_scores = F) + 
  scale_x_continuous(limits = c(-5000, 5000), breaks = seq(-4000, 4000, 2000), exp= c(0.005, 0.005), labels = scales::comma) +
  geom_hline(yintercept = 0, linetype = "dashed") + 
  labs(title = "2: Convert correlation values to Methodical scores")  
ggsave(plot = tubb6_methodical_scores_plot, "~/promoter_project/presentation_figures/methodical_part2.pdf", width = 16, height = 9)

# Add smoothed curve to TUBB6 plot
tubb6_smoothed_methodical_scores_plot = plotMethodicalScores(meth_site_values = tubb6_cpg_correlations, 
  reference_tss = tubb6_tss, p_value_threshold = NULL, smooth_scores = T, smoothed_curve_colour = "hotpink2", curve_alpha = 1, xlabel = "Distance to TSS (bp)") +
geom_hline(yintercept = 0, linetype = "dashed") + 
  theme(legend.position = c(0.9, 0.15)) + guides(fill = "none") +
  labs(title = "3: Smooth Methodical scores using exponential moving average") +
  scale_x_continuous(limits = c(-5000, 5000), breaks = seq(-4000, 4000, 2000), exp= c(0.005, 0.005), labels = scales::comma) 
ggsave(plot = tubb6_smoothed_methodical_scores_plot, "~/promoter_project/presentation_figures/methodical_part3.pdf", width = 16, height = 9)

# Add TMRs and significance thresholds to plot
tubb6_tmr_plot = plotTMRs(meth_site_plot = tubb6_smoothed_methodical_scores_plot, tmrs_gr = tubb6_tmrs, reference_tss = tubb6_tss) +
  geom_hline(yintercept = 0, linetype = "dashed") + 
  labs(title = "4: Identify TMRs where smoothed Methodical curve crosses significance thresholds") +
  geom_hline(yintercept = log10(0.005), linetype = "dashed", colour = "#7B5C90") +
  geom_hline(yintercept = -log10(0.005), linetype = "dashed", colour = "#BFAB25") +
  guides(fill = "none")
tubb6_tmr_plot2 = plotTMRs(meth_site_plot = tubb6_smoothed_methodical_scores_plot, tmrs_gr = tubb6_tmrs, reference_tss = tubb6_tss) +
  geom_hline(yintercept = 0, linetype = "dashed") + 
  labs(title = "4: Identify clusters where smoothed Methodical curve crosses significance thresholds") +
  geom_hline(yintercept = log10(0.005), linetype = "dashed", colour = "#7B5C90") +
  geom_hline(yintercept = -log10(0.005), linetype = "dashed", colour = "#BFAB25") +
  guides(fill = "none") + labs(color = "Cluster Direction")
ggsave(plot = tubb6_tmr_plot2, "~/promoter_project/presentation_figures/methodical_part4.pdf", width = 16, height = 9)

# Combine the plots 
tubb6_workflow_list = list(tubb6_cpg_correlation_plot, tubb6_methodical_scores_plot, tubb6_smoothed_methodical_scores_plot, tubb6_tmr_plot)
tubb6_workflow_list = lapply(tubb6_workflow_list, function(x) x + theme(plot.title = element_text(hjust = 0.5, size = 16)))
tubb6_workflow_plot = ggarrange(plotlist = tubb6_workflow_list, nrow = 2, ncol = 2, align = "hv")
ggsave(plot = tubb6_workflow_plot, "tmr_example_plots/methodical_workflow_plot.pdf", width = 20.57, height = 13.5)

# Create a function which will plot the CpG correlations for a specified transcript
plot_cpg_values_for_transcript = function(transcript, column_name = "cor", ylabel = "DNA Methylation-Transcription Correlation"){
  
  plot = plotMethSiteCorCoefs(meth_site_cor_values = dplyr::rename(cpgea_normal_correlation_results[[transcript]], meth_site = "cpg_name"),
    reference_tss = plyranges::filter(tss_gr, transcript_id == transcript), value_colours = "set2", 
    xlabel = sprintf("Distance to *%s* TSS (bp)", plyranges::filter(tss_gr, transcript_id == transcript)$gene_name),  
    ylabel = ylabel, title = "") +
    geom_hline(yintercept = 0, linetype = "dashed") +
    scale_x_continuous(limits = c(-5000, 5000), breaks = seq(-4000, 4000, 2000), exp= c(0.005, 0.005), labels = scales::comma) +
    theme(axis.title.x = ggtext::element_markdown(hjust = 0.5, size = 20))
  
  return(plot)
  
}

# Show correlation plots with TMRs for OTX1 and RHOF
otx1_tmr_plot = plotTMRs(meth_site_plot = plot_cpg_values_for_transcript("ENST00000366671"), 
  tmrs_gr = cpgea_normal_tmrs, reference_tss = tss_gr, transcript_id = "ENST00000366671") + 
  theme(legend.position = c(0.9, 0.15)) + guides(fill = "none")
rhof_tmr_plot = plotTMRs(meth_site_plot = plot_cpg_values_for_transcript("ENST00000267205"), 
  tmrs_gr = cpgea_normal_tmrs, reference_tss = tss_gr, transcript_id = "ENST00000267205") +
  theme(legend.position = c(0.9, 0.85)) + guides(fill = "none")
rhof_tmr_plot2 = plotTMRs(meth_site_plot = plot_cpg_values_for_transcript("ENST00000267205", ylabel = "Methylation-Transcription Correlation"), 
  tmrs_gr = cpgea_normal_tmrs, reference_tss = tss_gr, transcript_id = "ENST00000267205") +
  theme(legend.position = c(0.9, 0.85)) + guides(fill = "none") + labs(color = "Cluster Direction")
pitx2_tmr_plot = plotTMRs(meth_site_plot = plot_cpg_values_for_transcript("ENST00000394595"), 
  tmrs_gr = cpgea_normal_tmrs, reference_tss = tss_gr, transcript_id = "ENST00000394595") +
  theme(legend.position = c(0.9, 0.85)) + guides(fill = "none")
cracr2b_tmr_plot = plotTMRs(meth_site_plot = plot_cpg_values_for_transcript("ENST00000525077"), 
  tmrs_gr = cpgea_normal_tmrs, reference_tss = tss_gr, transcript_id = "ENST00000525077") +
  theme(legend.position = c(0.925, 0.85)) + guides(fill = "none")
cracr2b_tmr_plot2 = plotTMRs(meth_site_plot = plot_cpg_values_for_transcript("ENST00000525077", ylabel = "Methylation-Transcription Correlation"), 
  tmrs_gr = cpgea_normal_tmrs, reference_tss = tss_gr, transcript_id = "ENST00000525077") +
  theme(legend.position = c(0.925, 0.85)) + guides(fill = "none") + labs(color = "Cluster Direction")

# Get correlation values for clusters
cpgea_normal_clusters_correlations = data.table::fread("tmr_evaluation/cpgea_normal_tmrs_5kb_correlations_normal_samples.tsv.gz")
cpgea_normal_clusters_correlations$significance = plotR::sig_sym(cpgea_normal_clusters_correlations$q_val, symbol = "\u204E")

# Create a function which will plot the correlations for TMRs associated with a specified transcript
plot_tmr_correlations = function(transcript, title = "TMR Correlations"){
  
  # Get TMRs for transcript
  tmrs = plyranges::filter(cpgea_normal_tmrs, transcript_id == transcript)
  
  # Get strand of TSS associated with TMRs
  tss_strand = as.character(strand(GRanges(tmrs$tss_location)[1]))
  
  # Sort TMRs
  tmrs = sort(tmrs, decreasing = ifelse(tss_strand == "-", T, F))
  
  # Add correlation
  tmrs$correlation = cpgea_normal_clusters_correlations$cor[match(tmrs$tmr_name, cpgea_normal_clusters_correlations$genomic_region_name)]
  tmrs$significance = cpgea_normal_clusters_correlations$significance[match(tmrs$tmr_name, cpgea_normal_clusters_correlations$genomic_region_name)]
  tmrs_df = data.frame(tmrs)
  
  # Rename TMRs
  tmrs_df$tmr_name = paste("TMR", seq_along(tmrs))
  
  # Create plot
  ggplot(tmrs_df, 
    aes(y = correlation, x = tmr_name, fill = direction, label = significance)) + 
    geom_col(color = "black") + 
    geom_text(vjust = ifelse(tmrs_df$correlation >= 0, "bottom", "top"), size = 8) + 
    geom_hline(yintercept = 0) +
    scale_fill_manual(values = plotR::colour_list$purple_and_gold_light) +
    theme_classic() +
    scale_y_continuous(exp= expansion(mult = c(0.05, 0.05))) +
    labs(x = NULL, y = title, 
      title = NULL) +
    theme(axis.title.y = ggtext::element_markdown(hjust = 0.5, size = 18), 
      axis.text = element_text(size = 14), legend.position = "None", 
      ) 

}

# Create barplots with TMRs correlation 
rhof_tmr_correlation_barplot = plot_tmr_correlations(transcript = "ENST00000267205", title = "*RHOF* TMR Correlations")
rhof_tmr_correlation_barplot2 = plot_tmr_correlations(transcript = "ENST00000267205", title = "*RHOF* Cluster Correlations") +
  scale_x_discrete(labels = paste("Cluster", 1:3, sep = "\n")) + theme(axis.text = element_text(size = 11))
cracr2b_tmr_correlation_barplot = plot_tmr_correlations("ENST00000525077", title = "*CRACR2B* TMR Correlations") 
cracr2b_tmr_correlation_barplot2 = plot_tmr_correlations("ENST00000525077", title = "*CRACR2B* Cluster Correlations")  + 
  scale_x_discrete(labels = paste("Cluster", 1:3, sep = "\n")) + theme(axis.text = element_text(size = 11))

# Create presentation plots for RHOF and CRACR2B
rhof_combined_plot = ggarrange(plotlist = list(rhof_tmr_plot2, rhof_tmr_correlation_barplot2), 
  ncol = 2, nrow = 1, align = "hv", widths = c(3.5, 1))
cracr2b_combined_plot = ggarrange(plotlist = list(cracr2b_tmr_plot2, cracr2b_tmr_correlation_barplot2), 
  ncol = 2, nrow = 1, align = "hv", widths = c(3.5, 1))
ggsave(plot = rhof_combined_plot, filename = "~/promoter_project/presentation_figures/rhof_tmrs_plot.pdf", device = cairo_pdf, width = 16, height = 9)
ggsave(plot = cracr2b_combined_plot, filename = "~/promoter_project/presentation_figures/cracr2b_tmrs_plot.pdf", device = cairo_pdf, width = 16, height = 9)

# Combine the RHOF and CRACR2B TMR plots to show as examples, one as a standalone plot and another to be combined with tubb6_workflow_plot
tmr_example_plots_standalone = ggarrange(plotlist = list(rhof_tmr_plot2, rhof_tmr_correlation_barplot, cracr2b_tmr_plot, cracr2b_tmr_correlation_barplot), 
  ncol = 2, nrow = 2, align = "hv", widths = c(3.5, 1, 3.5, 1), labels = c("A", "", "B", ""))
ggsave(plot = tmr_example_plots_standalone, "tmr_example_plots/tmr_example_plots.pdf", width = 20.57, height = 18, device = cairo_pdf)
tmr_example_plots = ggarrange(plotlist = list(rhof_tmr_plot, rhof_tmr_correlation_barplot, cracr2b_tmr_plot, cracr2b_tmr_correlation_barplot), 
  ncol = 2, nrow = 2, align = "hv", widths = c(3.5, 1, 3.5, 1), labels = c("B", "", "C", ""))

# Load TMR distribution plot
tmr_distribution_plot = readRDS("tmr_distribution_plots/tmrs_5kb_bins_plot_final.rds")

# Combine tubb6_workflow_plot and tmr_plots and save
final_workflow_plot = ggarrange(plotlist = list(tubb6_workflow_plot, tmr_example_plots, tmr_distribution_plot), nrow = 3, heights = c(1.5, 2, 0.75), labels = c("A", "", "D"))
ggsave(final_workflow_plot, filename = "tmr_example_plots/figure3.1.pdf", width = 20.57, height = 38.25, device = cairo_pdf)

### Presentation plots
# Create a theme for use with poster
poster_theme = theme(plot.title = element_text(hjust = 0.5, size = 18), 
	axis.title = element_text(size = 18), axis.text = element_text(size = 18), strip.text.x = element_text(size = 18),
	legend.text = element_text(size = 18), legend.title = element_text(size = 18), legend.key.size = unit(1, "cm")) 

# Combine the plots 
tubb6_workflow_list_poster = lapply(tubb6_workflow_list, function(x) x + poster_theme)
tubb6_workflow_plot_poster = ggarrange(plotlist = tubb6_workflow_list_poster, nrow = 2, ncol = 2, align = "hv")
ggsave(plot = tubb6_workflow_plot_poster, "~/promoter_project/martin_plots/methodical_workflow_plot.pdf", width = 20.57, height = 13.5)

### Make plots for tumours

# Load CPGEA tumour correlation results
cpgea_tumour_correlation_results = readRDS("~/promoter_project/meth_transcript_correlations/meth_transcript_cor_lists/transcript_meth_cors_cpgea_tumour_samples_5kb.rds")

# Gt GRanges for TMRs for tumour prostate samples
cpgea_tumour_tmrs = readRDS("~/promoter_project/final_tmrs/final_tmr_granges/cpgea_tumour_tmrs_5kb.rds")
cpgea_tumour_tmrs = readRDS("~/promoter_project/final_tmrs/tmr_granges/cpgea_normal_tmrs_5kb.rds")

# Create a function which will plot the CpG correlations for a specified transcript
plot_cpg_values_for_transcript_tumour = function(transcript, column_name = "cor", ylabel = "DNA Methylation-Transcription Correlation"){
  
  plot = plotMethSiteCorCoefs(meth_site_cor_values = dplyr::rename(cpgea_tumour_correlation_results[[transcript]], meth_site = "cpg_name"),
    reference_tss = plyranges::filter(tss_gr, transcript_id == transcript), value_colours = "set2", 
    xlabel = sprintf("Distance to *%s* TSS (bp)", plyranges::filter(tss_gr, transcript_id == transcript)$gene_name),  
    ylabel = ylabel, title = "") +
    geom_hline(yintercept = 0, linetype = "dashed") +
    scale_x_continuous(limits = c(-5000, 5000), breaks = seq(-4000, 4000, 2000), exp= c(0.005, 0.005), labels = scales::comma) +
    theme(axis.title.x = ggtext::element_markdown(hjust = 0.5, size = 20))
  
  return(plot)
  
}

ENST00000343478; 
trans = sample(cpgea_tumour_tmrs$transcript_id, 1)
plotTMRs(meth_site_plot = plot_cpg_values_for_transcript_tumour(trans), 
  tmrs_gr = cpgea_tumour_tmrs, reference_tss = tss_gr, transcript_id = trans) +
  theme(legend.position = c(0.9, 0.85)) + guides(fill = "none")

gstp1_tmrs = find_TMRs(dplyr::rename(cpgea_tumour_correlation_results$ENST00000398606, meth_site = "cpg_name"))

gstp1_methodical_scores_plot = plotMethodicalScores(
  meth_site_values = dplyr::rename(cpgea_tumour_correlation_results$ENST00000398606, meth_site = "cpg_name"), reference_tss = tss_gr[tss_gr$transcript_id == "ENST00000398606"], 
  xlabel = "Distance to TSS (bp)", smooth_scores = T, smoothed_curve_colour = "hotpink2", curve_alpha = 1) + 
  scale_x_continuous(limits = c(-5000, 5000), breaks = seq(-4000, 4000, 2000), exp= c(0.005, 0.005), labels = scales::comma) +
  geom_hline(yintercept = 0, linetype = "dashed") 
gstp1_methodical_scores_plot = methodical::plotTMRs(gstp1_methodical_scores_plot, tmrs_gr = gstp1_tmrs, reference_tss = tss_gr[tss_gr$transcript_id == "ENST00000398606"]) + 
  theme(legend.position = c(0.9, 0.85), legend.text = element_text(size = 16), 
    legend.title = element_text(size = 16), legend.key.size = unit(1, units = "cm")) +  
  guides(fill = "none") + labs(color = "Correlation Direction")
ggsave(plot = gstp1_methodical_scores_plot, filename = "gstp1_tmrs.pdf", width = 16, height = 9)

trans_list = "ENST00000398606"
plot_list = lapply(trans_list, function(trans)
  plotTMRs(meth_site_plot = plot_cpg_values_for_transcript_tumour(trans), 
  tmrs_gr = cpgea_tumour_tmrs, reference_tss = tss_gr, transcript_id = trans) +
  theme(legend.position = c(0.9, 0.85)) + guides(fill = "none"))
plotR::pdf_save(plotlist = trans_list, filename = "test.pdf")
