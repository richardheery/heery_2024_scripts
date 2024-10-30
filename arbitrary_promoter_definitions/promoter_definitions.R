# Create a list with different promoter definitions for protein-coding transcripts 

# Promoter definitions come from the following 5 papers:
# A: Chen = -200-0 bp. From Tissue-independent and tissue-specific patterns of DNA methylation alteration in cancer; Epigenetics Chromatin. 2016
# B: Saghafinia -300-300 bp. From Pan-Cancer Landscape of Aberrant DNA Methylation across Human Tumors; Cell Reports. 2018
# C: Li = -2000-200 bp. From A genomic and epigenomic atlas of prostate cancer in Asian populations; Nature 2020
# D: Noushmehr -1500-1500 bp. From Identification of a CpG Island Methylator Phenotype that Defines a Distinct Subgroup of Glioma; Cancer Cell. 2010
# E: Cao = -4500-500 bp. From Multi-faceted epigenetic dysregulation of gene expression promotes esophageal squamous cell carcinoma; Nature Communications. 2020

# Load required packages
library(ggplot2)

# Make a data.frame with various promoter definitions ordered from smallest to biggest
promoter_definition_df = 
  data.frame(
    upstream = c(200, 300, 2000, 1500, 4500),
    downstream = c(0,  300, 200, 1500, 500),
    definition = c("A", "B", "C", "D", "E")
    )
promoter_definition_df$definition = factor(promoter_definition_df$definition, levels = rev(promoter_definition_df$definition))

# Save table with promoter definitions
saveRDS(promoter_definition_df, "promoter_definition_df.rds")

# Get TSS sites for protein-coding transcripts
pcg_transcript_tss_range = readRDS("../auxillary_data/genomic_annotation/pc_tss_gr.rds")

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