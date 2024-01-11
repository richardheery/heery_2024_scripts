# Calculate correlations between promoter methylation transcript expression for different promoter definitions.

# Load required packages
library(dplyr)
library(ggplot2)
library(plotR)
library(ggpubr)
library(doParallel)
library(cowplot)
library(methodical)

# Get names of MANE transcripts
mane_transcripts = gsub("\\.[0-9]*", "", readRDS("~/genomes/gencode/gencode_granges/gencode_v38_mane_transcript_ids.rds"))

# Get paths to all promoter methylation definition tables
promoter_definition_methylation_tables = readRDS("promoter_definition_methylation_tables/promoter_definition_methylation_tables.rds")

# Get kallisto output for protein-coding genes subset for normal tumour samples 
kallisto_deseq2_normalized_counts_pcg = data.frame(data.table::fread("~/mounts/local_mount/wgbs/cpgea/rnaseq/kallisto_tables/kallisto_deseq2_normalized_counts_pcg.tsv.gz"), row.names = 1)
kallisto_deseq2_normalized_counts_pcg_normal = dplyr::select(kallisto_deseq2_normalized_counts_pcg, starts_with("N"))
kallisto_deseq2_normalized_counts_pcg_tumour = dplyr::select(kallisto_deseq2_normalized_counts_pcg, starts_with("T"))

# Calculate correlation values between promoter methylation transcript expression for different promoter definitions in normal samples. Took 50 minutes. 
system.time({for(definition in names(promoter_definition_methylation_tables)){
  methylation_table = promoter_definition_methylation_tables[[definition]]
  feature_matches_df = data.frame(cluster = row.names(methylation_table), row.names(methylation_table))
  correlation_results = correlateR:::cor_tables(table1 = methylation_table, table2 = kallisto_deseq2_normalized_counts_pcg_normal, 
    feature_matches = feature_matches_df, calc_significance = T)
  data.table::fwrite(correlation_results, 
    paste0("promoter_definition_transcript_correlation_tables/", definition, "_definition_normal_sample_correlations.tsv.gz"), 
    sep = "\t", row.names = F, quote = F)
}})

# Calculate correlation values between promoter methylation transcript expression for different promoter definitions in tumour samples. Took 50 minutes. 
system.time({for(definition in names(promoter_definition_methylation_tables)){
  methylation_table = promoter_definition_methylation_tables[[definition]]
  feature_matches_df = data.frame(cluster = row.names(methylation_table), row.names(methylation_table))
  correlation_results = correlateR:::cor_tables(table1 = methylation_table, table2 = kallisto_deseq2_normalized_counts_pcg_tumour, 
    feature_matches = feature_matches_df, calc_significance = T)
  data.table::fwrite(correlation_results, 
    paste0("promoter_definition_transcript_correlation_tables/", definition, "_definition_tumour_sample_correlations.tsv.gz"), 
    sep = "\t", row.names = F, quote = F)
}})

### Make plots for normal samples

# Get lists of all normal correlation tables
normal_sample_correlation_tables = list.files("promoter_definition_transcript_correlation_tables", pattern = "normal_sample_correlations", full.names = T)
names(normal_sample_correlation_tables) = LETTERS[1:5]

# Create a list with all normal sample correlation tables
normal_sample_correlation_list = lapply(normal_sample_correlation_tables, function(x) 
  setNames(data.table::fread(x, select = 1:5), c("table1_feature", "table2_feature", "cor", "p_value", "q_value")))

# Combine the list into a single table
normal_sample_correlation_tables_combined = bind_rows(normal_sample_correlation_list, .id = "definition")
normal_sample_correlation_tables_combined$definition = 
  factor(normal_sample_correlation_tables_combined$definition, levels = LETTERS[1:5])

# Remove correlations where q_value is NA
normal_sample_correlation_tables_combined = filter(normal_sample_correlation_tables_combined, !is.na(q_value))

# Convert q-value to a significance symbol
normal_sample_correlation_tables_combined$significance = plotR::sig_sym(normal_sample_correlation_tables_combined$q_value, symbol = "\u204E")

# Make violin plots for distributions of correlation values without clusters
normal_sample_all_correlations_violin_plots = 
  ggplot(normal_sample_correlation_tables_combined, aes(y = cor, x =  definition, fill = definition)) +
    geom_violin() +
    geom_hline(yintercept = 0, linetype = "dashed") +
    theme_bw() +
    theme(plot.title = element_text(hjust = 0.5, size = 24), axis.title = element_text(size = 20), 
    axis.text = element_text(size = 14), strip.text.x = element_text(size = 14), legend.position = "None") +
    scale_fill_manual(values = c(colour_list$nine_greens[c(2, 4, 6, 8, 9)])) + 
    scale_y_continuous(limits = c(-1, 1), exp= c(0, 0), breaks = seq(-1, 1, 0.25)) +
    labs(x = "Promoter Definition", y = "Spearman Correlation", 
      title =  NULL)

# Denote whether correlations are positive, negative or uncorrelated
normal_sample_correlation_tables_combined = mutate(normal_sample_correlation_tables_combined, 
  correlation = case_when(
    q_value < 0.05 & cor > 0 ~ "Positive",
    q_value < 0.05 & cor < 0 ~ "Negative",
    q_value > 0.05 ~ "Uncorrelated"
    )
  )

# Convert correlation to a factor
normal_sample_correlation_tables_combined$correlation = factor(normal_sample_correlation_tables_combined$correlation, levels = c("Negative", "Uncorrelated", "Positive"))

# Get proportion of hypermethylated, hypomethylated unchanged promoters for each definition
normal_sample_correlation_tables_combined_summary = mutate(
  summarize(group_by(normal_sample_correlation_tables_combined, definition, correlation), count = n()),
  freq = count/sum(count))

# Create a barplot of proportion of hypermethylated, hypomethylated unchanged promoters for each definition
normal_correlation_proportions_barplot = ggplot(filter(normal_sample_correlation_tables_combined_summary, correlation != "Uncorrelated"), aes(y = freq, x = definition, fill = correlation)) +
 geom_col(position = "dodge", color  = "black")

# Adjust theme of barplot and save
normal_correlation_proportions_barplot = customize_ggplot_theme(normal_correlation_proportions_barplot, 
  title = "Proportion of Statistically Significant Correlations\nin Normal Prostate", 
  xlab = "Promoter Definition", ylab = "Proportion of Significant Correlations", fill_colors = c("#7B5C90B4", "#bfab25B4"), 
  plot_title_size = 28, axis_title_size = 24, axis_text_size = 20, legend_title_size = 24, legend_text_size = 20, legend_key_size = 1.5) + 
  scale_y_continuous(limits = c(0, 0.25), exp = c(0, 0))
ggsave(plot = normal_correlation_proportions_barplot, filename = "../presentation_figures/normal_correlation_proportions_barplot.pdf", width = 16, height = 9)

### Repeat for tumour samples

# Get lists of all tumour correlation tables
tumour_sample_correlation_tables = list.files("promoter_definition_transcript_correlation_tables", pattern = "tumour_sample_correlations", full.names = T)
names(tumour_sample_correlation_tables) = c("E", "A", "C", "D", "B")

# Create a list with all tumour sample correlation tables
tumour_sample_correlation_list = lapply(tumour_sample_correlation_tables, function(x) 
  setNames(data.table::fread(x, select = 1:5), c("table1_feature", "table2_feature", "cor", "p_value", "q_value")))

# Combine the list into a single table
tumour_sample_correlation_tables_combined = bind_rows(tumour_sample_correlation_list, .id = "definition")
tumour_sample_correlation_tables_combined$definition = 
  factor(tumour_sample_correlation_tables_combined$definition, levels = LETTERS[1:5])

# Remove correlations where q_value is NA
tumour_sample_correlation_tables_combined = filter(tumour_sample_correlation_tables_combined, !is.na(q_value))

# Make violin plots for distributions of correlation values without clusters
tumour_sample_all_correlations_violin_plots = 
  ggplot(tumour_sample_correlation_tables_combined, aes(y = cor, x =  definition, fill = definition)) +
    geom_violin() +
    geom_hline(yintercept = 0, linetype = "dashed") +
    theme_bw() +
    theme(plot.title = element_text(hjust = 0.5, size = 24), axis.title = element_text(size = 20), 
    axis.text = element_text(size = 14), strip.text.x = element_text(size = 14), legend.position = "None") +
    scale_fill_manual(values = c(colour_list$nine_greens[c(2, 4, 6, 8, 9)])) + 
    scale_y_continuous(limits = c(-1, 1), exp= c(0, 0), breaks = seq(-1, 1, 0.25)) +
    labs(x = "Promoter Definition", y = "Spearman Correlation", 
      title = NULL)

# Denote whether correlations are positive, negative or uncorrelated
tumour_sample_correlation_tables_combined = mutate(tumour_sample_correlation_tables_combined, 
  correlation = case_when(
    q_value < 0.05 & cor > 0 ~ "Positive",
    q_value < 0.05 & cor < 0 ~ "Negative",
    q_value > 0.05 ~ "Uncorrelated"
    )
  )

# Convert correlation to a factor
tumour_sample_correlation_tables_combined$correlation = factor(tumour_sample_correlation_tables_combined$correlation, levels = c("Negative", "Uncorrelated", "Positive"))

# Get proportion of hypermethylated, hypomethylated unchanged promoters for each definition
tumour_sample_correlation_tables_combined_summary = mutate(
  summarize(group_by(tumour_sample_correlation_tables_combined, definition, correlation), count = n()),
  freq = count/sum(count))

# Create a barplot of proportion of hypermethylated, hypomethylated unchanged promoters for each definition
tumour_correlation_proportions_barplot = ggplot(filter(tumour_sample_correlation_tables_combined_summary, correlation != "Uncorrelated"), aes(y = freq, x = definition, fill = correlation)) +
 geom_col(position = "dodge", color  = "black")

# Adjust theme of barplot save
tumour_correlation_proportions_barplot = customize_ggplot_theme(tumour_correlation_proportions_barplot, title = NULL, 
  xlab = "Promoter Definition", ylab = "Proportion of Significant Correlations", fill_colors = colour_list$purple_and_gold_light[c(1, 2)]) +
  scale_y_continuous(limits = c(0, 0.25), exp= c(0, 0)) + theme(legend.position = c(0.85, 0.9))

combined_tumour_plots = cowplot::plot_grid(plotlist = list(tumour_sample_all_correlations_violin_plots, 
  tumour_correlation_proportions_barplot), align = "h", nrow = 1, labels = c("A", "B"))
ggsave(plot = combined_tumour_plots, "promoter_definition_plots/correlation_plots/combined_tumour_plots.pdf", width = 18, height = 9)

combined_normal_plots = cowplot::plot_grid(plotlist = list(normal_sample_all_correlations_violin_plots, 
  normal_correlation_proportions_barplot + labs(title = NULL) + theme(legend.position = c(0.85, 0.85))), align = "h", nrow = 1, labels = c("D", "E"))

# Combine normal tumour correlation proportion tables
combined_sample_correlation_tables = bind_rows(
  Normal = normal_sample_correlation_tables_combined,
  Tumour = tumour_sample_correlation_tables_combined, .id = "group"
)

# Combine normal tumour correlation proportion violin_plots save
combined_correlation_violin_plots = ggplot(combined_sample_correlation_tables, aes(y = cor, x=  definition, fill = definition)) +
  geom_violin() +
  geom_hline(yintercept = 0, linetype = "dashed") +
  scale_y_continuous(limits = c(-1, 1), exp= c(0, 0), breaks = seq(-1, 1, 0.25)) 
combined_correlation_violin_plots = customize_ggplot_theme(combined_correlation_violin_plots, title = NULL, 
  xlab = "Promoter Definition", ylab = "DNA Methylation-Transcription Correlation", fill_colors = colour_list$nine_greens[c(2, 4, 6, 8, 9)], 
  facet = "group", facet_scales = "fixed", facet_labels = c("Normal Prostate", "Prostate Cancer"), show_legend = T) + theme(strip.background = element_blank())
combined_correlation_violin_plots = combined_correlation_violin_plots + 
  guides(fill = guide_legend(override.aes = list(alpha= 0 , color = "white"))) +
  theme(legend.text = element_text(color = "transparent"))

# Save one plot for paper and another for presentation 
ggsave(plot = combined_correlation_violin_plots + labs(title = NULL), "promoter_definition_plots/correlation_plots/combined_correlation_violin_plots.pdf", width = 16, height = 9, bg = "white")
ggsave(plot = combined_correlation_violin_plots + 
    labs(title = "Correlation Between Promoter Methylation and Transcript Expression"), "~/promoter_project/presentation_figures/combined_correlation_violin_plots.pdf", width = 16, height = 9, bg = "white")

# Combine normal tumour correlation proportion tables
combined_sample_correlation_tables_combined_summary = bind_rows(
  Normal = normal_sample_correlation_tables_combined_summary,
  Tumour = tumour_sample_correlation_tables_combined_summary, .id = "group"
)

# Combine normal tumour correlation proportion barplots save
combined_correlation_proportions_barplot = ggplot(filter(combined_sample_correlation_tables_combined_summary, correlation != "Uncorrelated"), aes(y = freq, x = definition, fill = correlation)) +
 geom_col(position = "dodge", color  = "black")
combined_correlation_proportions_barplot = customize_ggplot_theme(combined_correlation_proportions_barplot, title = "Proportion of Significant Correlations", 
  xlab = "Promoter Definition", ylab = "Proportion of Significant Correlations", fill_colors = c("#7B5C90B4", "#bfab25B4"), 
  facet = "group", facet_scales = "fixed", facet_labels = c("Normal Prostate", "Prostate Cancer")) + theme(strip.background = element_blank()) +
  scale_y_continuous(limits = c(0, 0.25), exp= c(0, 0))

# Save one plot for paper and another for presentation 
ggsave(plot = combined_correlation_proportions_barplot + labs(title = NULL), 
  "promoter_definition_plots/correlation_plots/combined_correlation_proportions_barplot_new.pdf", width = 16, height = 9)
ggsave(plot = combined_correlation_proportions_barplot + labs(y = "Proportion"), 
  "~/promoter_project/presentation_figures/combined_correlation_proportions_barplot.pdf", width = 16, height = 9)

# Extract the legend from barplot
barplot_legend = as_ggplot(get_legend(combined_correlation_proportions_barplot))

###

# Get GRanges for TSS sites
tss_gr = readRDS("~/genomes/gencode/gencode_granges/pcg_transcript_tss_ranges_gencode_v38.rds")

# Load promoter region plot
promoter_region_plot = readRDS("promoter_definition_plots/promoter_definition_plot.rds")

# Load CPGEA normal correlation results
cpgea_normal_correlation_results = readRDS("../meth_transcript_correlations/meth_transcript_cor_lists/transcript_meth_cors_cpgea_normal_samples_5kb.rds")

# Create a function which will plot the CpG correlations for a specified transcript
plot_cpg_values_for_transcript = function(transcript, column_name = "correlation", title = NULL, 
  ylabel = "DNA Methylation-Transcription Correlation", xlabel = "Distance of CpG to TSS (bp)"){
  
  plot = methodical::plotMethSiteCorCoefs(meth_site_cor_values = dplyr::rename(cpgea_normal_correlation_results[[transcript]], "meth_site" = "cpg_name"),
    reference_tss = plyranges::filter(tss_gr, transcript_id == transcript), value_colours = "set2", xlabel = xlabel,  
    ylabel = ylabel, title = title) +
    geom_hline(yintercept = 0, linetype = "dashed") +
    scale_x_continuous(limits = c(-5000, 5000), breaks = seq(-4000, 4000, 2000), exp= c(0.005, 0.005), labels = scales::comma) +
    theme(axis.title.x = ggtext::element_markdown(hjust = 0.5, size = 20))
  
  return(plot)
  
}

# Create a function which will plot the promoter correlations for a specified transcript
plot_promoter_correlations = function(transcript_id, title = "Promoter\nMethylation Change"){
  
  # Get differential methylation results for transcript
  transcript_results = dplyr::filter(normal_sample_correlation_tables_combined, table1_feature == transcript_id)
  
  # Create plot
  ggplot(transcript_results, 
    aes(y = cor, x = definition, fill = definition, label = significance)) + 
    geom_col(color = "black") + 
    geom_text(vjust = ifelse(transcript_results$cor >= 0, "bottom", "top"), size = 8) + 
    geom_hline(yintercept = 0) +
    scale_fill_manual(values = RColorBrewer::brewer.pal(9, "YlGn")[c(2, 4, 6, 8, 9)]) +
    theme_classic() +
    scale_y_continuous(exp= expansion(mult = c(0.05, 0.05))) +
    labs(x = "Promoter Definition", y = title, 
      title = NULL) +
    theme(plot.title = ggtext::element_markdown(hjust = 0.5, size = 18), 
      axis.title.y = ggtext::element_markdown(size = 20), axis.title.x = ggtext::element_markdown(size = 20),
      axis.text = element_text(size = 18), axis.text.x = element_text(hjust = 1), legend.position = "None")

}

# Create CpG and promoter correlation plots for RHOF
rhof_plot_cpg_correlations = plot_cpg_values_for_transcript("ENST00000267205", 
  xlabel = "Distance to *RHOF* TSS (bp)")
rhof_plot_cpg_correlations_pres = plot_cpg_values_for_transcript("ENST00000267205", 
  title = expression("Correlation of DNA Methylation with Transcription Near" ~ italic("RHOF") ~ "TSS in Normal Prostate"), ylabel = "Methylation-Transcription Correlation")
rhof_plot_promoter_correlations = plot_promoter_correlations(transcript_id = "ENST00000267205",
  title = "*RHOF* Promoter Methylation-<br>Transcription Correlation")

# Create CpG and promoter correlation plots for OTX1
otx1_plot_cpg_correlations = plot_cpg_values_for_transcript("ENST00000366671")
otx1_plot_promoter_correlations = plot_promoter_correlations(transcript_id = "ENST00000366671",
  title = "*OTX1* Promoter Methylation-<br>Transcription Correlation")

# Create CpG and promoter correlation plots for PITX1
pitx1_plot_cpg_correlations = plot_cpg_values_for_transcript("ENST00000394595")
pitx1_plot_promoter_correlations = plot_promoter_correlations(transcript_id = "ENST00000394595",
  title = "*PITX1* Promoter Methylation-<br>Transcription Correlation")

# # Create CpG and promoter correlation plots for TUBB6
tubb6_plot_cpg_correlations = plot_cpg_values_for_transcript("ENST00000591909")
tubb6_plot_promoter_correlations = plot_promoter_correlations(transcript_id = "ENST00000591909", 
  title = "*TUBB6* Promoter Methylation<br>Transcription Correlation")

# Create CpG and promoter correlation plots for CRACR2B
cracr2b_plot_cpg_correlations = plot_cpg_values_for_transcript("ENST00000525077", 
  xlabel = "Distance to *CRACR2B* TSS (bp)")
cracr2b_plot_promoter_correlations = plot_promoter_correlations(transcript_id = "ENST00000525077",
  title = "*CRACR2* Promoter Methylation-<br>Transcription Correlation")

# Create a plot just for CRACR2B and save
cracr2b_plotlist = list(cracr2b_plot_cpg_correlations, cracr2b_plot_promoter_correlations, 
  promoter_region_plot, NULL) 
cracr2b_combined_plot = plot_grid(plotlist = cracr2b_plotlist, nrow = 2, ncol = 2, align = "hv", 
  rel_heights = c(5.5, 2), rel_widths = c(3.5, 1))
ggsave(plot = cracr2b_combined_plot, filename = "~/promoter_project/presentation_figures/cracr2b_combined_meth_plot.pdf", 
  width = 20.57, height = 12.27, device = cairo_pdf)

promoter_correlation_plotlist = list(
  promoter_region_plot, NULL,
  rhof_plot_cpg_correlations, rhof_plot_promoter_correlations,
  cracr2b_plot_cpg_correlations, cracr2b_plot_promoter_correlations
  )

violin_and_barplot_list = list(
  combined_correlation_violin_plots,  
  combined_correlation_proportions_barplot + labs(title = NULL)
  )

violin_and_barplot = 
  plot_grid(plotlist = violin_and_barplot_list, nrow = 2, ncol = 1, align = "hv", 
    rel_widths = c(3.5, 1), rel_heights = c(5.5, 5.5), labels = c("D", "E"))

combined_correlation_plots = 
  plot_grid(plotlist = promoter_correlation_plotlist, nrow = 3, ncol = 2, align = "hv", 
    rel_widths = c(3.5, 1), rel_heights = c(1.5, 5.5, 5.5), labels = c("A", "", "B", "", "C", ""))
ggsave(plot = combined_correlation_plots, "promoter_definition_plots/correlation_plots/cpg_vs_promoter_correlation_plots2.pdf", 
  width = 20.57, height = 21.27, device = cairo_pdf)

# Combine all plots together
cpg_vs_promoter_correlation_panel_plot = plot_grid(combined_correlation_plots, violin_and_barplot, nrow = 2, rel_heights = c(12.5, 11))
cpg_vs_promoter_correlation_panel_plot = plot_grid(combined_correlation_plots, combined_normal_plots, nrow = 2, rel_heights = c(12.5, 5.5))

ggsave(plot = cpg_vs_promoter_correlation_panel_plot, filename = "promoter_definition_plots/correlation_plots/figure2.2.pdf", 
  width = 20.57, height = 30.08, device = cairo_pdf)

# Create a plot just for RHOF and save
rhof_plotlist = list(rhof_plot_cpg_correlations_pres, rhof_plot_promoter_correlations, 
  promoter_region_plot, NULL) 
rhof_combined_plot = plot_grid(plotlist = rhof_plotlist, nrow = 2, ncol = 2, align = "hv", 
  rel_heights = c(5.5, 2), rel_widths = c(3.5, 1))
ggsave(plot = rhof_combined_plot, filename = "~/promoter_project/presentation_figures/rhof_combined_correlation_plot.pdf", 
  width = 20.57, height = 12.27, device = cairo_pdf)

# Create a plot just for OTX1 and save
otx1_plotlist = list(otx1_plot_cpg_correlations, otx1_plot_promoter_correlations, 
  promoter_region_plot, NULL) 
otx1_combined_plot = plot_grid(plotlist = otx1_plotlist, nrow = 2, ncol = 2, align = "hv", 
  rel_heights = c(5.5, 2), rel_widths = c(3.5, 1))
ggsave(plot = otx1_combined_plot, filename = "~/promoter_project/presentation_figures/otx1_combined_correlation_plot.pdf", 
  width = 20.57, height = 12.27, device = cairo_pdf)

# Create a plot just for PITX2 and save
pitx2_plotlist = list(pitx2_plot_cpg_correlations, pitx2_plot_promoter_correlations, 
  promoter_region_plot, NULL) 
pitx2_combined_plot = plot_grid(plotlist = pitx2_plotlist, nrow = 2, ncol = 2, align = "hv", 
  rel_heights = c(5.5, 2), rel_widths = c(3.5, 1))
ggsave(plot = pitx2_combined_plot, filename = "~/promoter_project/presentation_figures/pitx2_combined_correlation_plot.pdf", 
  width = 20.57, height = 12.27, device = cairo_pdf)

### Presentation plots
normal_sample_all_correlations_violin_plots = 
  ggplot(normal_sample_correlation_tables_combined, aes(y = cor, x =  definition, fill = definition)) +
    geom_violin() +
    geom_hline(yintercept = 0, linetype = "dashed") +
    theme_bw() +
    theme(plot.title = element_text(hjust = 0.5, size = 28), axis.title = element_text(size = 24), 
    axis.text = element_text(size = 20), strip.text.x = element_text(size = 14), legend.position = "None") +
    scale_fill_manual(values = c(colour_list$nine_greens[c(2, 4, 6, 8, 9)])) + 
    scale_y_continuous(limits = c(-1, 1), exp= c(0, 0), breaks = seq(-1, 1, 0.25)) +
    labs(x = "Promoter Definition", y = "Correlation Value", 
      title =  "Correlation Between Promoter Methylation and\nTranscript Expression in Normal Prostate")
normal_sample_all_correlations_violin_plots
ggsave(plot = normal_sample_all_correlations_violin_plots, "../presentation_figures/normal_sample_all_correlations_violin_plots.pdf", width = 16, height = 9)

### Make plots for MANE transcripts

# Filter combined_sample_correlation_tables for MANE transcripts
combined_sample_correlation_tables_mane = filter(combined_sample_correlation_tables, table1_feature %in% mane_transcripts)

# Combine normal tumour correlation proportion violin_plots save
combined_correlation_violin_plots_mane = ggplot(combined_sample_correlation_tables_mane, aes(y = cor, x=  definition, fill = definition)) +
  geom_violin() +
  geom_hline(yintercept = 0, linetype = "dashed") +
  scale_y_continuous(limits = c(-1, 1), exp= c(0, 0), breaks = seq(-1, 1, 0.25)) 
combined_correlation_violin_plots_mane = customize_ggplot_theme(combined_correlation_violin_plots_mane, title = NULL, 
  xlab = "Promoter Definition", ylab = "DNA Methylation-Transcription Correlation", fill_colors = colour_list$nine_greens[c(2, 4, 6, 8, 9)], 
  facet = "group", facet_scales = "fixed", facet_labels = c("Normal Prostate", "Prostate Cancer"), show_legend = T) + theme(strip.background = element_blank())
combined_correlation_violin_plots_mane = combined_correlation_violin_plots_mane + 
  guides(fill = guide_legend(override.aes = list(alpha= 0 , color = "white"))) +
  theme(legend.text = element_text(color = "transparent"))
ggsave(plot = combined_correlation_violin_plots_mane, 
  "promoter_definition_plots/mane_plots/combined_correlation_violin_plots.pdf", width = 16, height = 9)

# Get proportion of hypermethylated, hypomethylated unchanged promoters for each definition
normal_sample_correlation_tables_combined_summary_mane = mutate(
  summarize(group_by(filter(normal_sample_correlation_tables_combined, table1_feature %in% mane_transcripts), definition, correlation), count = n()),
  freq = count/sum(count))
tumour_sample_correlation_tables_combined_summary_mane = mutate(
  summarize(group_by(filter(tumour_sample_correlation_tables_combined, table1_feature %in% mane_transcripts), definition, correlation), count = n()),
  freq = count/sum(count))

# Combine normal tumour correlation proportion tables
combined_sample_correlation_tables_combined_summary_mane = bind_rows(
  Normal = normal_sample_correlation_tables_combined_summary_mane,
  Tumour = tumour_sample_correlation_tables_combined_summary_mane, .id = "group"
)

# Combine normal tumour correlation proportion barplots save
combined_correlation_proportions_barplot_mane = ggplot(filter(combined_sample_correlation_tables_combined_summary_mane, correlation != "Uncorrelated"), aes(y = freq, x = definition, fill = correlation)) +
 geom_col(position = "dodge", color  = "black")
combined_correlation_proportions_barplot_mane = customize_ggplot_theme(combined_correlation_proportions_barplot_mane, title = "Proportion of Significant Correlations", 
  xlab = "Promoter Definition", ylab = "Proportion of Significant Correlations", fill_colors = c("#7B5C90B4", "#bfab25B4"), 
  facet = "group", facet_scales = "fixed", facet_labels = c("Normal Prostate", "Prostate Cancer")) + theme(strip.background = element_blank()) +
  scale_y_continuous(limits = c(0, 0.25), exp= c(0, 0))
ggsave(plot = combined_correlation_proportions_barplot_mane, 
  "promoter_definition_plots/mane_plots/combined_correlation_proportions_barplot.pdf", width = 16, height = 9)
