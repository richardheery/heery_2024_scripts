# Make plots showing the impact of promoter definition on differential promoter methylation

# Load required packaged
library(dplyr)
library(methodical)
library(matrixStats)
library(cowplot)
library(ggpubr)

# Get TSS for all protein-coding transcripts
gencode_tss_gr = readRDS("~/genomes/gencode/gencode_granges/pcg_transcript_tss_ranges_gencode_v38.rds")

# Get methrix object for CPGEA samples
cpgea_meth_rse = HDF5Array::loadHDF5SummarizedExperiment("/media/rich/Elements/wgbs/cpgea/wgbs/cpgea_wgbs_hg38/")

# Create a data.frame with the coordinates of the different promoter definitions
promoter_definition_df = 
  data.frame(
    upstream = -c(200, 300, 2000, 1500, 4500),
    downstream = c(0,  300, 200, 1500, 500),
    definition = c("A", "B", "C", "D", "E")
    )
promoter_definition_df$definition = factor(promoter_definition_df$definition, levels = rev(promoter_definition_df$definition))

# Make promoter coordinates plot
promoter_region_plot = ggplot(promoter_definition_df, 
  aes(xmin = upstream, xmax = downstream, x = NULL, y = definition,  group = definition)) + 
  geom_linerange(linewidth = 10.5, position = position_dodge(0.06), color = RColorBrewer::brewer.pal(9, "YlGn")[c(2, 4, 6, 8, 9)]) +
  theme_classic() + 
  theme(plot.title = element_text(hjust = 0.5, size = 24), legend.text = element_text(size = 18),
    axis.title = element_text(size = 20), axis.text = element_text(size = 18))  +
  scale_x_continuous(limits = c(-5000, 5000), breaks = seq(-4000, 4000, 2000), expand = c(0.005, 0.005), labels = scales::comma) +
  geom_vline(xintercept = 0, linetype = "dashed") +
  labs(x = NULL, y = "Promoter\nDefinition")

# Get a list with the promoter methylation change values for all promoter definitions
promoter_definition_methylation_change = readRDS("promoter_definition_methylation_change.rds")

promoter_definition_differential_methylation_results = readRDS("promoter_differential_methylation_results/promoter_diff_methylation_results.rds")

# Combine promoter differential methylation results into a single table
promoter_definition_differential_methylation_results = bind_rows(promoter_definition_differential_methylation_results, .id = "definition")

# Convert q-value to a significance symbol
promoter_definition_differential_methylation_results$significance = plotR::sig_sym(promoter_definition_differential_methylation_results$q_value, symbol = "\u204E")

# Find promoters which have methylation values with all definitions
common_transcripts = names(which(table(unlist(lapply(promoter_definition_methylation_change, names))) == 5))

# Create a data.frame with the methylation values for all common_transcripts
promoter_definition_methylation_df = data.frame(lapply(promoter_definition_methylation_change, function(x) x[common_transcripts]))

# Find the minimum and maximum promoter methylation values for each transcript
min_methylation = rowMins(as.matrix(promoter_definition_methylation_df) , na.rm = T)
max_methylation = rowMaxs(as.matrix(promoter_definition_methylation_df) , na.rm = T)

# Identify transcripts which have methylation loss greater than 0.15 with some promoter definitions and methylation gain greater than 0.15 with others
high_variance_promoters = promoter_definition_methylation_df[which(min_methylation < -0.15 & max_methylation > 0.15), ]

# Filter for only those cases with methylation values for all promoter definitions
high_variance_promoters = high_variance_promoters[complete.cases(high_variance_promoters), ]

# Convert promoter_definition_methylation_df to long format
promoter_definition_methylation_df_long = reshape2::melt(tibble::rownames_to_column(promoter_definition_methylation_df, "transcript"))

# Create a function which plots promoter methylation change using the different promoter definitions
plot_promoter_methylation_change = function(transcript_id, title = "Promoter\nMethylation Change"){
  
  # Get differential methylation results for transcript
  transcript_results = dplyr::filter(promoter_definition_differential_methylation_results, feature == transcript_id)
  
  # Create plot
  ggplot(transcript_results, 
    aes(y = mean_change, x = definition, fill = definition, label = significance)) + 
    geom_col(color = "black") + 
    geom_text(vjust = ifelse(transcript_results$mean_change >= 0, "bottom", "top"), size = 8) + 
    geom_hline(yintercept = 0) +
    scale_fill_manual(values = RColorBrewer::brewer.pal(9, "YlGn")[c(2, 4, 6, 8, 9)]) +
    theme_classic() +
    scale_y_continuous(expand = expansion(mult = c(0.05, 0.05))) +
    labs(x = "Promoter Definition", y = title, 
      title = NULL) +
    theme(plot.title = ggtext::element_markdown(hjust = 0.5, size = 20), 
      axis.title.y = ggtext::element_markdown(size = 20), axis.title.x = element_text(size = 20),
      axis.text = element_text(size = 18), axis.text.x = element_text(hjust = 1), legend.position = "None") 

}

# Make a function which will plot CpG methylation change for all CpGs in +/- 5KB of a TSS
plot_cpg_methylation_change = function(transcript, title = NULL, xlabel = "Distance of CpG to TSS (bp)"){
  
  # Get TSS of transcript
  tss_gr = plyranges::filter(gencode_tss_gr, transcript_id == transcript)
  
  # Get methylation values of all CpG values within +/- 5 KB of the TSS
  transcript_methylation = methodical::extractGRangesMethSiteValues(cpgea_meth_rse, promoters(tss_gr, 5000, 5001))
  
  # Get mean methylation change for each site
  transcript_methylation_change = data.frame(meth_change = 
      rowMeans(select(transcript_methylation, starts_with("T")) 
        - select(transcript_methylation, starts_with("N")), na.rm = T))
  
  # Plot methylation change for the CpGs
  cpg_plot = methodical::plotMethylationValues(meth_site_values = transcript_methylation_change, 
    sample_name = "meth_change", reference_tss = tss_gr, 
    xlabel = xlabel, ylabel = "DNA Methylation Change", title = title) +
    geom_hline(yintercept = 0, linetype = "dashed") +
    scale_x_continuous(limits = c(-5000, 5000), breaks = seq(-4000, 4000, 2000), expand = c(0.005, 0.005), labels = scales::comma) +
    scale_y_continuous(breaks = seq(-1, 1, 0.2), expand = c(0.05, 0.05), labels = scales::comma)
  
}

# Make plots of CpG methylation change and promoter definition differential methylation for GNAL and TRO
gnal_cpg_meth_change_plot = plot_cpg_methylation_change(transcript = "ENST00000590228", 
  xlabel = expression("Distance to" ~ italic("GNAL") ~ "TSS (bp)"))
tro_cpg_meth_change_plot = plot_cpg_methylation_change(transcript = "ENST00000452830", 
  xlabel = expression("Distance to" ~ italic("TRO") ~ "TSS (bp)"))
gnal_promoters_plot = plot_promoter_methylation_change("ENST00000590228", 
  title = "*GNAL* Promoter<br>Methylation Change")
tro_promoters_plot = plot_promoter_methylation_change("ENST00000452830", 
  title = "*TRO* Promoter<br>Methylation Change")

# Combine plots into a single figure along with promoter definition plot and save
combined_plotlist = list(promoter_region_plot, NULL, 
  gnal_cpg_meth_change_plot, gnal_promoters_plot, 
  tro_cpg_meth_change_plot, tro_promoters_plot)
combined_plot = plot_grid(plotlist = combined_plotlist, nrow = 3, ncol = 2, align = "hv", 
  rel_heights = c(1.5, 5.5, 5.5), rel_widths = c(3.5, 1), labels = c("A", "", "B", "", "C", ""))
ggsave(plot = combined_plot, filename = "promoter_definition_plots/differential_methylation_plots/promoter_definition_meth_change_examples.pdf", 
  width = 20.57, height = 21.27, device = cairo_pdf)

# Create a plot just for GNAL and save
gnal_plotlist = list(gnal_cpg_meth_change_plot, gnal_promoters_plot, 
  promoter_region_plot, NULL) 
gnal_combined_plot = plot_grid(plotlist = gnal_plotlist, nrow = 2, ncol = 2, align = "hv", 
  rel_heights = c(5.5, 2), rel_widths = c(3.5, 1)) 
ggsave(plot = gnal_combined_plot, filename = "~/promoter_project/presentation_figures/gnal_combined_meth_plot.pdf", 
  width = 20.57, height = 12.27, device = cairo_pdf)

# Create a plot just for TRO and save
tro_plotlist = list(tro_cpg_meth_change_plot, tro_promoters_plot, 
  promoter_region_plot, NULL) 
tro_combined_plot = plot_grid(plotlist = tro_plotlist, nrow = 2, ncol = 2, align = "hv", 
  rel_heights = c(5.5, 2), rel_widths = c(3.5, 1))
ggsave(plot = tro_combined_plot, filename = "~/promoter_project/presentation_figures/tro_combined_meth_plot.pdf", 
  width = 20.57, height = 12.27, device = cairo_pdf)

# Load barplot with proportion of differentally methylated promoters for each definition
meth_change_proportions_barplot = readRDS("promoter_definition_plots/differential_methylation_plots/meth_change_proportions_barplot.rds")

# Add label "D" to barplot
meth_change_proportions_barplot = plot_grid(plotlist = list(meth_change_proportions_barplot), labels = "D")

# Add barplot to complete plot
complete_plotlist = list(combined_plot, meth_change_proportions_barplot)
complete_plot = plot_grid(plotlist = complete_plotlist, nrow = 2, ncol = 1,
  rel_heights = c(0.65, 0.35), labels = c("", "D"))
ggsave(plot = complete_plot, filename = "promoter_definition_plots/differential_methylation_plots/figure1.1.pdf", 
  width = 20.57, height = 32.72, device = cairo_pdf)
