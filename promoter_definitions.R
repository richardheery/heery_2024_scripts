# Create a list with different promoter definitions for all protein-coding transcripts annotated by Gencode and 
# create a plot showing the locations of these definitions relative to the TSS

# Definitions come from the following 5 papers:
# A: Chen = -200-0 bp. From Tissue-independent and tissue-specific patterns of DNA methylation alteration in cancer; Epigenetics Chromatin. 2016
# B: Saghafinia -300-300 bp. From Pan-Cancer Landscape of Aberrant DNA Methylation across Human Tumors; Cell Reports. 2018
# C: Li = -2000-200 bp. From A genomic and epigenomic atlas of prostate cancer in Asian populations; Nature 2020
# D: Noushmehr -1500-1500 bp. From Identification of a CpG Island Methylator Phenotype that Defines a Distinct Subgroup of Glioma; Cancer Cell. 2010
# E: Cao = -4500-500 bp. From Multi-faceted epigenetic dysregulation of gene expression promotes esophageal squamous cell carcinoma; Nature Communications. 2020

# Load required packages
library(methodicalFinal)
library(ggplot2)

# Make a data.frame with various promoter definitions ordered from smallest to biggest
promoter_definition_df = 
  data.frame(
    upstream = c(200, 300, 2000, 1500, 4500),
    downstream = c(0,  300, 200, 1500, 500),
    definition = c("A", "B", "C", "D", "E")
    )
promoter_definition_df$definition = factor(promoter_definition_df$definition, levels = rev(promoter_definition_df$definition))

# Get TSS sites for protein-coding transcripts
pcg_transcript_tss_range = readRDS("~/genomes/gencode/gencode_granges/pcg_transcript_tss_ranges_gencode_v38.rds")

# Create a list of promoters using the promoter definitions in promoter_definition_df
# Seems that TSS shold be considered 1
promoter_definition_list = lapply(1:nrow(promoter_definition_df), function(x)
  promoters(
    pcg_transcript_tss_range, 
    upstream = promoter_definition_df$upstream[x], 
    downstream = promoter_definition_df$downstream[x] + 1) 
)

# Set names for list
names(promoter_definition_list) = promoter_definition_df$definition

# Save list
saveRDS(promoter_definition_list, "promoter_definition_list.rds")

# Create a plot showing different promoter definitions
promoter_region_plot = ggplot(promoter_definition_df , aes(xmin = -upstream, xmax = downstream, x = NULL, y = definition,  group = definition)) + 
    geom_linerange(linewidth = 11, position = position_dodge(0.06), color = plotR::colour_list$nine_greens[c(2, 4, 6, 8, 9)]) +
    theme_classic() + 
    theme(plot.title = element_text(hjust = 0.5, size = 24), legend.text = element_text(size = 18),
      axis.title = element_text(size = 20), axis.text = element_text(size = 18))  +
    scale_x_continuous(limits = c(-5000, 5000), expand = c(0, 0), breaks = seq(-4000, 4000, 2000), labels = scales::comma) +
    labs(x = "Distance to TSS (bp)", y = "Promoter\nDefinition", title = NULL) +
  geom_vline(xintercept = 0, linetype = "dashed")
ggsave(plot = promoter_region_plot, "promoter_definition_plots/promoter_definition_plot.pdf", height = 2.25, width = 16)
ggsave(plot = promoter_region_plot + ggtitle("Location of Published Promoter Definitions Relative to TSS"), 
  "../presentation_figures/promoter_definition_plot.pdf", height = 2.75, width = 16)
saveRDS(promoter_region_plot, "promoter_definition_plots/promoter_definition_plot.rds")
