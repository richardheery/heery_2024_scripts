# Annotate TMRs using genomic features and chromatin states

# Load required packages
library(GenomicRanges)
library(dplyr)
library(patchwork)
source("../auxillary_scripts/plotting_functions.R")

# Get list of TMRs
tmr_list = readRDS("../finding_tmrs/tmr_list.rds")

# Combine tmr_list into a single GRanges
tmrs_combined = unlist(GRangesList(tmr_list))
tmrs_combined$group = names(tmrs_combined)

# Get a GRanges for transcripts and expand
transcripts_gr = readRDS("../auxillary_data/genomic_annotation/pc_transcripts_gr.rds")
transcripts_gr = methodical:::.expand_transcript_ranges(transcripts_gr, 5000, 5000)

# Get CpG islands from UCSC
cpg_island_annotation <- annotatr::build_annotations(genome = "hg38", annotations = "hg38_cpgs")
cpg_island_annotation$type = factor(gsub("hg38_", "", cpg_island_annotation$type))
cpg_island_annotation$type = recode_factor(cpg_island_annotation$type, 
  "cpg_inter" = "Open Sea",
  "cpg_islands" = "CpG Islands",
  "cpg_shelves" = "CpG Shelves",
  "cpg_shores" = "CpG Shores")
cpg_island_annotation$region_type = cpg_island_annotation$type
cpg_island_annotation$type = NULL

# Load hg38 genome annotation
genome_annotation_hg38 = readRDS("../auxillary_data/genomic_annotation/genome_annotation_hg38.rds")

# Remove CpG Islands
genome_annotation_hg38 = filter_gr(genome_annotation_hg38, "region_type", values = "CpG Island", invert = T)

# Shorten names of exons and introns
genome_annotation_hg38$region_type[genome_annotation_hg38$region_type == "exon"] = "Exons"
genome_annotation_hg38$region_type[genome_annotation_hg38$region_type == "intron"] = "Introns"

# Remove " Region" from cluster_annotation$region_type
genome_annotation_hg38$region_type = gsub(" Region", "", genome_annotation_hg38$region_type)

# Remove repeats from annotation
genome_annotation_hg38 = plyranges::filter(genome_annotation_hg38, !region_type %in% 
    c("Alu", "CR1", "DNA Transposon", "L1", "L2", "LTR", "Low Complexity", "MIR", "SVA", "Satellite", "Simple Repeat", "Protein-Coding"))

# Remove lncRNA regions
genome_annotation_hg38 = plyranges::filter(genome_annotation_hg38, !grepl("lncRNA", region_type))

# Combine genome_annotation_hg38 and cpg_island_annotation
genome_annotation_hg38 = c(genome_annotation_hg38, cpg_island_annotation)

# Add transcripts_gr as background TMR search space
background_gr = transcripts_gr
mcols(background_gr) = NULL
background_gr$region_type = "Background"
genome_annotation_hg38 = c(genome_annotation_hg38, background_gr)

# Convert genome_annotation_hg38 to a list
genome_annotation_hg38_list = GRangesList(split(genome_annotation_hg38, genome_annotation_hg38$region_type))

# Calculate the size of the overlaps of transcripts_gr with each chromatin state
annotation_widths = sapply(genome_annotation_hg38_list, function(x) calculate_regions_intersections(gr1 = transcripts_gr, gr2 = x))

# Find overlaps between TMRs and chromatin states
tmr_annotation_overlaps = data.frame(findOverlaps(tmrs_combined, genome_annotation_hg38))
tmr_annotation_overlaps$tmr_group = gsub("_tmrs_.*", "_tmrs", tmrs_combined$group[tmr_annotation_overlaps$queryHits])
tmr_annotation_overlaps$direction = tmrs_combined$direction[tmr_annotation_overlaps$queryHits]
tmr_annotation_overlaps$annotation = genome_annotation_hg38$region_type[tmr_annotation_overlaps$subjectHits]

# Make a table summarizing the overlaps
tmr_annotation_overlaps_summary = dplyr::summarise(group_by(tmr_annotation_overlaps, tmr_group, direction, annotation), count = n())

# Normalize count by MB covered by the chromatin state
tmr_annotation_overlaps_summary$normalized_count = tmr_annotation_overlaps_summary$count/annotation_widths[tmr_annotation_overlaps_summary$annotation]*1e6

# Convert region_type to a factor and give specified order
tmr_annotation_overlaps_summary = filter(tmr_annotation_overlaps_summary, annotation %in% 
    c("Background", "CpG Islands", "CpG Shores", "CpG Shelves", "Open Sea", "Predicted Promoter", "Predicted Enhancer", "Open Chromatin", "CTCF BS"))
tmr_annotation_overlaps_summary$annotation = factor(tmr_annotation_overlaps_summary$annotation, 
  c("Background", "CpG Islands", "CpG Shores", "CpG Shelves", "Open Sea", "Predicted Promoter", "Predicted Enhancer", "Open Chromatin", "CTCF BS"))

# Add data set name to tmr_annotation_overlaps_summary
tmr_annotation_overlaps_summary$tmr_group = with(tmr_annotation_overlaps_summary, case_when(
  tmr_group == "cpgea_normal_tmrs" ~ "Normal Prostate",
  tmr_group == "cpgea_tumour_tmrs" ~ "Prostate Tumours",
  tmr_group == "mcrpc_tmrs" ~ "Prostate Metastases"
))

# Put tmr_groups in right order
tmr_annotation_overlaps_summary$tmr_group = factor(tmr_annotation_overlaps_summary$tmr_group, 
  levels = c("Normal Prostate", "Prostate Tumours", "Prostate Metastases"))

# Create a plot annotating TMRs and save
cluster_genomic_feature_annotation_plot = ggplot(tmr_annotation_overlaps_summary, aes(x = tmr_group, y = normalized_count, fill = direction)) +
  geom_col(position = "dodge", colour = "black") 
cluster_genomic_feature_annotation_plot = customize_ggplot_theme(cluster_genomic_feature_annotation_plot, title = NULL, 
  xlab = "Dataset", ylab = "Number of TMRs per MB", fill_title = "TMR Direction", fill_colors = colour_list$purple_and_gold_light, x_labels_angle = 55, colors = "black",
  facet = "annotation", facet_nrow = 1, facet_scales = "free_x", strip_text_size = 20, axis_text_size = 14, 
  legend_title_size = 24, legend_text_size = 20, legend_key_size = 1.5) + 
  theme(strip.background = element_blank(), plot.margin = margin(t = 0.5, unit = "cm")) +
  guides(linetype = guide_legend(override.aes = list(legend.text = element_text(size = 50))))
cluster_genomic_feature_annotation_plot

### Annotate chromatin states

# Get a GRanges with chromatin states for prostate 
prostate_18_states_hg38_gr = readRDS("../auxillary_data/genomic_annotation/prostate_18_states_hg38_gr.rds")

# Shorten "PolyComb" in description to "PC" 
levels(prostate_18_states_hg38_gr$description)[3] = "Flanking TSS 5'"
levels(prostate_18_states_hg38_gr$description)[4] = "Flanking TSS 3'"
levels(prostate_18_states_hg38_gr$description)[5] = "Strong Transcription"
levels(prostate_18_states_hg38_gr$description)[6] = "Weak Transcription"
levels(prostate_18_states_hg38_gr$description)[7] = "Genic Enhancer 1"
levels(prostate_18_states_hg38_gr$description)[8] = "Genic Enhancer 2"
levels(prostate_18_states_hg38_gr$description)[12] = "ZNF Genes/Repeats"
levels(prostate_18_states_hg38_gr$description)[16] = "Repressed PC"
levels(prostate_18_states_hg38_gr$description)[17] = "Weak Repressed PC"

# Split prostate_18_states_hg38_gr into a list
prostate_18_states_hg38_gr$region_type = prostate_18_states_hg38_gr$description
background_gr2 = background_gr
background_gr2$region_type = "Background2"
prostate_18_states_hg38_gr = c(background_gr, prostate_18_states_hg38_gr, background_gr2)
prostate_18_states_hg38_gr$region_type = factor(prostate_18_states_hg38_gr$region_type, 
  levels = c(levels(prostate_18_states_hg38_gr$region_type)[1:10], levels(prostate_18_states_hg38_gr$region_type)[20], 
    levels(prostate_18_states_hg38_gr$region_type)[11:19]))
prostate_18_states_hg38_gr_list = split(prostate_18_states_hg38_gr, prostate_18_states_hg38_gr$region_type)

# Calculate the size of the ovelaps of transcripts_gr with each chromatin state
chromatin_state_widths = sapply(prostate_18_states_hg38_gr_list, function(x) calculate_regions_intersections(gr1 = transcripts_gr, gr2 = x))

# Find overlaps between TMRs and chromatin states
tmr_state_overlaps = data.frame(findOverlaps(tmrs_combined, prostate_18_states_hg38_gr))
tmr_state_overlaps$tmr_group = gsub("_tmrs_.*", "_tmrs", tmrs_combined$group[tmr_state_overlaps$queryHits])
tmr_state_overlaps$direction = tmrs_combined$direction[tmr_state_overlaps$queryHits]
tmr_state_overlaps$chromatin_state = prostate_18_states_hg38_gr$region_type[tmr_state_overlaps$subjectHits]

# Make a table summarizing the overlaps
tmr_state_overlaps_summary = summarise(group_by(tmr_state_overlaps, tmr_group, direction, chromatin_state), count = n())

# Normalize count by MB covered by the chromatin state
tmr_state_overlaps_summary$normalized_count = tmr_state_overlaps_summary$count/chromatin_state_widths[tmr_state_overlaps_summary$chromatin_state]*1e6

# Add data set name to tmr_state_overlaps_summary
tmr_state_overlaps_summary$tmr_group = with(tmr_state_overlaps_summary, case_when(
  tmr_group == "cpgea_normal_tmrs" ~ "Normal Prostate",
  tmr_group == "cpgea_tumour_tmrs" ~ "Prostate Tumours",
  tmr_group == "mcrpc_tmrs" ~ "Prostate Metastases"
))

# Put tmr_group in right order
tmr_state_overlaps_summary$tmr_group = factor(tmr_state_overlaps_summary$tmr_group, 
  levels = c("Normal Prostate", "Prostate Tumours", "Prostate Metastases"))

# Create a plot annotating TMRs and save
cluster_chromatin_state_annotation_plot = ggplot(tmr_state_overlaps_summary, aes(x = tmr_group, y = normalized_count, fill = direction)) +
  geom_col(position = "dodge", colour = "black")
cluster_chromatin_state_annotation_plot = customize_ggplot_theme(cluster_chromatin_state_annotation_plot, title = NULL, 
  xlab = "Dataset", ylab = "Number of TMRs per MB", fill_title = "TMR Direction", fill_colors = colour_list$purple_and_gold_light, x_labels_angle = 55, colors = "black",
  facet = "chromatin_state", facet_nrow = 2, facet_scales = "free_x", strip_text_size = 17, axis_text_size = 16, axis_title_size = 24, legend_key_size = 1.5, 
  legend_title_size = 24, legend_text_size = 20) + 
  theme(strip.background = element_blank(), plot.margin = margin(t = 0.5, unit = "cm")) +
  guides(linetype = guide_legend(override.aes = list(legend.text = element_text(size = 50))))
cluster_chromatin_state_annotation_plot

### Combine The TMR distribution plots and TMR annotation plots 
combined_annotation_plot = cluster_genomic_feature_annotation_plot / cluster_chromatin_state_annotation_plot +
  plot_layout(heights = c(1, 2), nrow = 2, ncol = 1) +
  plot_annotation(tag_levels = 'A') & theme(plot.tag = element_text(size = 20))
ggsave(plot = combined_annotation_plot, filename = "../figures/figure5A_and_B.pdf", width = 27, height = 27)