# Calculate correlations between promoter methylation and transcript expression for different promoter definitions.

# Load required packages
library(dplyr)
library(ggplot2)
library(ggpubr)
library(doParallel)
library(cowplot)
library(methodical)

# Get paths to all promoter methylation definition tables. 
promoter_definition_methylation_tables = readRDS("promoter_definition_methylation_tables/promoter_definition_methylation_tables.rds")

# Get kallisto output for protein-coding genes subset for normal tumour samples 
kallisto_deseq2_normalized_counts_pcg = data.frame(data.table::fread("../auxillary_data/rnaseq_data/cpgea_kallisto_deseq2_normalized_counts_pcg.tsv.gz"), row.names = 1)
kallisto_deseq2_normalized_counts_pcg_normal = dplyr::select(kallisto_deseq2_normalized_counts_pcg, starts_with("N"))
kallisto_deseq2_normalized_counts_pcg_tumour = dplyr::select(kallisto_deseq2_normalized_counts_pcg, starts_with("T"))

# Create a directory for promoter methylation-transcription correlation results
dir.create("promoter_definition_transcript_correlation_tables")

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

# Calculate correlation values between promoter methylation transcript expression for different promoter definitions in metastases samples. Took 50 minutes. 
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
normal_sample_correlation_tables_combined$significance = ::sig_sym(normal_sample_correlation_tables_combined$q_value, symbol = "\u204E")

# Make violin plots for distributions of correlation values without clusters
normal_sample_all_correlations_violin_plots = 
  ggplot(normal_sample_correlation_tables_combined, aes(y = cor, x =  definition, fill = definition)) +
    geom_violin() +
    geom_hline(yintercept = 0, linetype = "dashed") +
    theme_bw() +
    theme(plot.title = element_text(hjust = 0.5, size = 28), axis.title = element_text(size = 24), 
    axis.text = element_text(size = 20), strip.text.x = element_text(size = 14), legend.position = "None") +
    scale_fill_manual(values = c(colour_list$nine_greens[c(2, 4, 6, 8, 9)])) + 
    scale_y_continuous(limits = c(-1, 1), exp= c(0, 0), breaks = seq(-1, 1, 0.25)) +
    labs(x = "Promoter Definition", y = "Promoter Methylation vs\nTranscription Correlation Values", 
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
  title = NULL, 
  xlab = "Promoter Definition", ylab = "Proportion of Significant Correlations", fill_colors = c("#7B5C90B4", "#bfab25B4"), 
  plot_title_size = 28, axis_title_size = 24, axis_text_size = 20, legend_title_size = 24, legend_text_size = 20, legend_key_size = 1.5) + 
  scale_y_continuous(limits = c(0, 0.1), exp = c(0, 0)) + theme(legend.position = c(0.85, 0.9))

# Combine normal plots and save
combined_normal_correlation_plots = ggarrange(plotlist = list(normal_sample_all_correlations_violin_plots, normal_correlation_proportions_barplot),
  labels = c("D", "E"))
ggsave(plot = combined_normal_correlation_plots, filename = "../figures/figureD_and_E.pdf", width = 20.57, height = 12.17)

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
    labs(x = "Promoter Definition", y = "Promoter Methylation vs\nTranscription Correlation Values", 
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
  xlab = "Promoter Definition", ylab = "Proportion of Significant Correlations", fill_colors = colour_list$purple_and_gold_light) +
  theme(legend.position = c(0.9, 0.9))

# Combine tumour plots and save
combined_tumour_plots = ggarrange(plotlist = list(tumour_sample_all_correlations_violin_plots, tumour_correlation_proportions_barplot),
  labels = c("A", "B"))
ggsave(plot = combined_tumour_plots, filename = "../final_plots/supp_figure7.pdf", width = 20.57, height = 12.17)

