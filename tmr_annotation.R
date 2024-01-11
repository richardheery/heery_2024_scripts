# Annotate TMRs using genomic features and chromatin states

# Load required packages
library(methodical)
library(dplyr)
library(plotR)
library(patchwork)

# Get list of TMRs
tmr_list = readRDS("final_tmr_granges/final_tmr_5kb_list.rds")

# Get a GRanges for TSS sites
tss_gr = readRDS("~/genomes/gencode/gencode_granges/pcg_transcript_tss_ranges_gencode_v38.rds")
# tss_gr_mane = tss_gr[tss_gr$transcript_id %in% gencode_transcript_annotation_mane$ID]
# tss_proximal_ranges_mane = promoters(tss_gr_mane, 0, 1001)

# Get ranges within +/- 5 KB of TSS
tss_proximal_ranges = promoters(tss_gr, 5000, 5001)

# Load hg38 genome annotation
genome_annotation_hg38 = methodicalFinal::genome_annotation_hg38

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

genome_annotation_hg38 = GRangesList(split(genome_annotation_hg38, genome_annotation_hg38$region_type))

# Annotate tss_proximal_ranges
background_genomic_feature_annotation = annotateGRanges(genomic_regions = tss_proximal_ranges, annotation_ranges = genome_annotation_hg38, overlap_measure = "proportion")

# Convert background_genomic_feature_annotation into a data.frame
background_genomic_feature_annotation = data.frame(region_type = names(background_genomic_feature_annotation),
  proportion = background_genomic_feature_annotation, background = "Background", row.names = NULL)

# Annotate clusters
cluster_genomic_feature_annotation = lapply(tmr_list, function(x) 
  annotateGRanges(genomic_regions = x, annotation_ranges = genome_annotation_hg38, ignore.strand = T, overlap_measure = "proportion"))

# Combine cluster annotation into a single data.frame and convert to long format
cluster_genomic_feature_annotation = bind_rows(cluster_genomic_feature_annotation, .id = "cluster_type")
cluster_genomic_feature_annotation = tidyr::pivot_longer(cluster_genomic_feature_annotation, -cluster_type, names_to = "region_type", values_to = "proportion")

# Get the max proportion overlap for each type of element and select only those where the max is at least 0.014 (to include open chromatin).
max_proportion = data.frame(summarise(group_by(cluster_genomic_feature_annotation, region_type), max_proportion = max(proportion)))
selected_region_types = filter(max_proportion, max_proportion >= 0.014)$region_type
cluster_genomic_feature_annotation_selected = filter(cluster_genomic_feature_annotation, region_type %in% selected_region_types)

# Convert region_type to a factor and give specified order
cluster_genomic_feature_annotation_selected$region_type = factor(cluster_genomic_feature_annotation_selected$region_type, 
  c("CpG Island", "Predicted Promoter", "Predicted Enhancer", "Open Chromatin", "CTCF BS", "Exons", "Introns"))

# Filter background_genomic_feature_annotation for regions in cluster_genomic_feature_annotation_selected
background_genomic_feature_annotation = dplyr::filter(background_genomic_feature_annotation, region_type %in% cluster_genomic_feature_annotation_selected$region_type)
background_genomic_feature_annotation$region_type = factor(background_genomic_feature_annotation$region_type, levels = levels(cluster_genomic_feature_annotation_selected$region_type))

# Add data set name to cluster_genomic_feature_annotation_selected
cluster_genomic_feature_annotation_selected$dataset = gsub("_negative|_positive", "", cluster_genomic_feature_annotation_selected$cluster_type)
cluster_genomic_feature_annotation_selected$dataset = with(cluster_genomic_feature_annotation_selected, case_when(
  dataset == "cpgea_normal_tmrs" ~ "Normal Prostate",
  dataset == "cpgea_tumour_tmrs" ~ "Prostate Tumours",
  dataset == "mcrpc_tmrs" ~ "Prostate Metastases"
))

# Put datasets in right order
cluster_genomic_feature_annotation_selected$dataset = factor(cluster_genomic_feature_annotation_selected$dataset, 
  levels = c("Normal Prostate", "Prostate Tumours", "Prostate Metastases"))

# Add direction to cluster_genomic_feature_annotation_selected
cluster_genomic_feature_annotation_selected$direction = stringr::str_to_title(gsub(".*_", "", cluster_genomic_feature_annotation_selected$cluster_type))

# Create a plot annotating TMRs and save
cluster_genomic_feature_annotation_plot = ggplot(filter(cluster_genomic_feature_annotation_selected, cluster_type != "background"), aes(x = dataset, y = proportion, fill = direction)) +
  geom_col(position = "dodge", colour = "black") + 
  geom_hline(data = background_genomic_feature_annotation, mapping = aes(yintercept = proportion, color = background), linetype = "dashed")
cluster_genomic_feature_annotation_plot = customize_ggplot_theme(cluster_genomic_feature_annotation_plot, title = NULL, 
  xlab = "Dataset", ylab = "Overlap with Genomic Features", fill_title = "TMR Direction", fill_colors = colour_list$purple_and_gold_light, x_labels_angle = 55, colors = "black",
  facet = "region_type", facet_nrow = 1, facet_scales = "free_x", strip_text_size = 20, axis_text_size = 14, 
  legend_title_size = 24, legend_text_size = 20, legend_key_size = 1.5) + 
  theme(strip.background = element_blank(), plot.margin = margin(t = 0.5, unit = "cm")) +
  guides(linetype = guide_legend(override.aes = list(legend.text = element_text(size = 50))))
cluster_genomic_feature_annotation_plot
ggsave(plot = cluster_genomic_feature_annotation_plot, filename = "tmr_annotation_plots/cluster_genomic_feature_annotation_plot.pdf", width = 16, height = 9)

# cluster_genomic_feature_annotation_selected$region_type = factor(cluster_genomic_feature_annotation_selected$region_type, 
#   levels = c("Predicted Promoter", "Exons", "Introns", "CpG Island", "Predicted Enhancer"))
# cluster_genomic_feature_annotation_selected = filter(cluster_genomic_feature_annotation_selected, !is.na(region_type))
# background_genomic_feature_annotation = filter(background_genomic_feature_annotation, region_type %in% cluster_genomic_feature_annotation_selected$region_type)
# cluster_genomic_feature_annotation_plot = ggplot(filter(cluster_genomic_feature_annotation_selected, cluster_type != "background"), aes(x = dataset, y = proportion, fill = direction)) +
#   geom_col(position = "dodge", colour = "black") + 
#   geom_hline(data = background_genomic_feature_annotation, mapping = aes(yintercept = proportion, color = background), linetype = "dashed")
# cluster_genomic_feature_annotation_plot = customize_ggplot_theme(cluster_genomic_feature_annotation_plot, title = NULL, 
#   xlab = "Dataset", ylab = "Overlap with Genomic Features", fill_title = "Cluster\nDirection", fill_colors = colour_list$purple_and_gold_light, x_labels_angle = 55, colors = "black",
#   facet = "region_type", facet_nrow = 1, facet_scales = "free_x", strip_text_size = 18, axis_text_size = 16) + 
#   theme(strip.background = element_blank(), plot.margin = margin(t = 0.5, unit = "cm")) +
#   guides(linetype = guide_legend(override.aes = list(legend.text = element_text(size = 50))))
# cluster_genomic_feature_annotation_plot
ggsave(plot = cluster_genomic_feature_annotation_plot, "../presentation_figures/cluster_genomic_feature_annotation_plot.pdf", width = 16, height = 9)

### Annotate chromatin states

# Get a GRanges with chromatin states for prostate
prostate_18_states_hg38_gr = readRDS("~/roadmap/prostate_roadmap/prostate_18_states_hg38_gr.rds")
prostate_18_states_hg38_gr = split(prostate_18_states_hg38_gr, prostate_18_states_hg38_gr$description)

# Shorten "PolyComb" in description to "PC" 
levels(prostate_18_states_hg38_gr$description)[3] = "Flanking TSS 5'"
levels(prostate_18_states_hg38_gr$description)[4] = "Flanking TSS 3'"
levels(prostate_18_states_hg38_gr$description)[3] = "Flanking TSS 5'"
levels(prostate_18_states_hg38_gr$description)[12] = "ZNF Genes/Repeats"
levels(prostate_18_states_hg38_gr$description)[16] = "Repressed PC"
levels(prostate_18_states_hg38_gr$description)[17] = "Weak Repressed PC"

# Annotate tss_proximal_ranges
background_chromatin_state_annotation = annotateGRanges(genomic_regions = tss_proximal_ranges, 
  annotation_ranges = prostate_18_states_hg38_gr, overlap_measure = "proportion")

# Convert background_chromatin_state_annotation into a data.frame
background_chromatin_state_annotation = data.frame(region_type = names(background_chromatin_state_annotation),
  proportion = background_chromatin_state_annotation, background = "Background", row.names = NULL)

# Put levels in order of prostate_18_states_hg38_gr
background_chromatin_state_annotation$region_type = factor(background_chromatin_state_annotation$region_type, 
  levels = names(prostate_18_states_hg38_gr))

# Annotate clusters
cluster_chromatin_state_annotation = lapply(tmr_list, function(x) 
  annotateGRanges(genomic_regions = x, prostate_18_states_hg38_gr, overlap_measure = "proportion"))

# Combine cluster chromatin_state_annotation into a single data.frame and convert to long format
cluster_chromatin_state_annotation = bind_rows(cluster_chromatin_state_annotation, .id = "cluster_type")
cluster_chromatin_state_annotation = tidyr::pivot_longer(cluster_chromatin_state_annotation, -cluster_type, names_to = "region_type", values_to = "proportion")

# Put levels in order of prostate_18_states_hg38_gr
cluster_chromatin_state_annotation$region_type = factor(cluster_chromatin_state_annotation$region_type, 
  levels = names(prostate_18_states_hg38_gr))

# Add data set name to cluster_chromatin_state_annotation
cluster_chromatin_state_annotation$dataset = gsub("_negative|_positive", "", cluster_chromatin_state_annotation$cluster_type)
cluster_chromatin_state_annotation$dataset = with(cluster_chromatin_state_annotation, case_when(
  dataset == "cpgea_normal_tmrs" ~ "Normal Prostate",
  dataset == "cpgea_tumour_tmrs" ~ "Prostate Tumours",
  dataset == "mcrpc_tmrs" ~ "Prostate Metastases"
))

# Put datasets in right order
cluster_chromatin_state_annotation$dataset = factor(cluster_chromatin_state_annotation$dataset, 
  levels = c("Normal Prostate", "Prostate Tumours", "Prostate Metastases"))

# Add direction to cluster_chromatin_state_annotation
cluster_chromatin_state_annotation$direction = stringr::str_to_title(gsub(".*_", "", cluster_chromatin_state_annotation$cluster_type))

# Create a plot annotating TMRs and save
cluster_chromatin_state_annotation_plot = ggplot(filter(cluster_chromatin_state_annotation, cluster_type != "background"), 
  aes(x = dataset, y = proportion, fill = direction)) +
  geom_col(position = "dodge", colour = "black") + 
  geom_hline(data = background_chromatin_state_annotation, mapping = aes(yintercept = proportion, color = background), linetype = "dashed")
cluster_chromatin_state_annotation_plot = customize_ggplot_theme(cluster_chromatin_state_annotation_plot, title = NULL, 
  xlab = "Dataset", ylab = "Overlap with Chromatin States", fill_title = "TMR Direction", fill_colors = colour_list$purple_and_gold_light, x_labels_angle = 55, colors = "black",
  facet = "region_type", facet_nrow = 2, facet_scales = "free_x", strip_text_size = 20, axis_text_size = 16, axis_title_size = 24, legend_key_size = 1.5, 
  legend_title_size = 24, legend_text_size = 20) + 
  theme(strip.background = element_blank(), plot.margin = margin(t = 0.5, unit = "cm")) +
  guides(linetype = guide_legend(override.aes = list(legend.text = element_text(size = 50))))
cluster_chromatin_state_annotation_plot
ggsave(plot = cluster_chromatin_state_annotation_plot, filename = "tmr_annotation_plots/cluster_chromatin_state_annotation_plot.pdf", width = 27*1, height = 16*1)

### Combine The TMR distribution plots and TMR annotation plots 

# # Load the TMR bins plot
# tmrs_5kb_bins_plot = readRDS("tmr_distribution_plots/tmrs_5kb_bins_plot_final.rds")

# Combine plots with patchwork
# patchwork_plot = tmrs_5kb_bins_plot / cluster_genomic_feature_annotation_plot / cluster_chromatin_state_annotation_plot +
#   plot_layout(heights = c(1, 1, 2), nrow = 3, ncol = 1) +
#   plot_annotation(tag_levels = 'A') & theme(plot.tag = element_text(size = 20))
# ggsave(plot = patchwork_plot, filename = "tmr_annotation_plots/figure4.pdf", width = 27, height = 36)

patchwork_plot = cluster_genomic_feature_annotation_plot / cluster_chromatin_state_annotation_plot +
  plot_layout(heights = c(1, 2), nrow = 2, ncol = 1) +
  plot_annotation(tag_levels = 'A') & theme(plot.tag = element_text(size = 20))
ggsave(plot = patchwork_plot, filename = "tmr_annotation_plots/figure4.1.pdf", width = 27, height = 27)

# Create a theme for use with poster
poster_theme = theme(plot.title = element_text(hjust = 0.5, size = 30), 
	axis.title = element_text(size = 26), axis.text = element_text(size = 24), strip.text.x = element_text(size = 26),
	legend.text = element_text(size = 24), legend.title = element_text(size = 26), legend.key.size = unit(2, "cm")) 

patchwork_plot_poster = (cluster_genomic_feature_annotation_plot + poster_theme) / 
  (cluster_chromatin_state_annotation_plot + poster_theme + facet_wrap("region_type", nrow = 3)) + 
  plot_layout(heights = c(1, 2), nrow = 2, ncol = 1) +
  plot_annotation(tag_levels = 'A') & theme(plot.tag = element_text(size = 20))
ggsave(plot = patchwork_plot_poster, filename = "~/promoter_project/presentation_figures/tmr_annotation_plot_poster.pdf", width = 27, height = 27)

# Presentation plots
ggsave(plot = cluster_genomic_feature_annotation_plot + poster_theme, "~/promoter_project/martin_plots/cluster_genomic_feature_annotation_plot.pdf", width = 27, height = 16)
ggsave(plot = cluster_chromatin_state_annotation_plot + poster_theme + facet_wrap("region_type", nrow = 3), "~/promoter_project/martin_plots/cluster_chromatin_state_annotation_plot.pdf", width = 27, height = 16)
