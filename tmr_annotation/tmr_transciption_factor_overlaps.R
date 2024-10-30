# Evaluate the overlaps of clusters with transcriptional regulator (TR) binding sites

# Load required packages
library(GenomicRanges)
library(methodical)
library(doParallel)
source("../auxillary_scripts/enrichment_tests.R")

# Get a GRanges for transcripts and expand
transcripts_gr = readRDS("../auxillary_data/genomic_annotation/pc_transcripts_gr.rds")
transcripts_gr = methodical:::.expand_transcript_ranges(transcripts_gr, 5000, 5000)

# Get all CpG sites as a GRanges object
cpg_sites = readRDS("~/genomes/cpgs/all_cpg_sites_hg38.rds")

# Get all TSS cpgs
tss_cpgs = subsetByOverlaps(cpg_sites, transcripts_gr)
rm(cpg_sites); gc()

# Get a list with all TMR types
tmr_list = readRDS("../finding_tmrs/tmr_list.rds")

# For each TMR group get the TSS regions for the corresponding transcripts
tmr_list_transcript_regions = lapply(tmr_list, function(x) transcripts_gr[transcripts_gr$ID %in% x$ID])

# Load a list of GRanges for TR binding sites for prostate from ReMap
remap_prostate_gr_list = readRDS("../auxillary_data/genomic_annotation/remap_prostate_gr_list.rds")
all_binding_sites_gr = unlist(GRangesList(remap_prostate_gr_list))

# Calculate enrichment of TMRs for all binding sites
system.time({all_tf_enrichment_results = lapply(tmr_list, function(x) query_subject_cpg_overlap_test(query_regions = x, 
  subject_regions = all_binding_sites_gr, control_regions = transcripts_gr, background_cpgs = tss_cpgs, alternative = "greater"))})
all_tf_enrichment_results = bind_rows(all_tf_enrichment_results, .id = "tmr_group")
write.table(all_tf_enrichment_results, "all_tf_enrichment_results.tsv", sep = "\t", quote = F)

# Calculate enrichment of TF binding sites for TMRs compared to non-TMR regions using binding sites in prostate from ReMap and save. 
system.time({tmr_remap_enrichment = 
  foreach(tmr = names(tmr_list), .packages = "GenomicRanges") %do% {
    print(tmr)
    query_subject_cpg_overlap_test_apply(
      query_regions = tmr_list[[tmr]], control_regions = tmr_list_transcript_regions[[tmr]], 
      background_cpgs = tss_cpgs, gr_list = remap_prostate_gr_list, alternative = "two.sided")
  }
names(tmr_remap_enrichment) = names(tmr_list)
})
saveRDS(tmr_remap_enrichment, "tmr_remap_enrichment.rds")
tmr_remap_enrichment = readRDS("tmr_remap_enrichment.rds")

# Update names of tmr_remap_enrichment
names(tmr_remap_enrichment) = c("Normal Prostate Negative", "Normal Prostate Positive", "Prostate Tumours Negative", 
  "Prostate Tumours Positive", "Prostate Metastases Negative", "Prostate Metastases Positive")

# Get all significantly enriched results
tmr_remap_enrichment_significant = lapply(tmr_remap_enrichment, function(x) 
  dplyr::filter(x, q_value < 0.05, enrichment > 1))

# Save all significant results as a single table
helpR::fwrite(bind_rows(tmr_remap_enrichment_significant, .id = "group"), 
  filename = "tf_bs_results.tsv")

# Get all significant TRs
significant_tfs = lapply(tmr_remap_enrichment_significant, function(x) x$tf_name)

# Get top 20 results for each TMR group
tmr_remap_enrichment_top = lapply(tmr_remap_enrichment_significant, function(x) head(arrange(x, desc(enrichment)), 20))

# Combine results for all TMRs
combined_tmr_remap_enrichment = bind_rows(tmr_remap_enrichment_top, .id = "group")

# Put group in desired order
combined_tmr_remap_enrichment$group = factor(combined_tmr_remap_enrichment$group, unique(combined_tmr_remap_enrichment$group))

# Update column names of combined_tmr_remap_enrichment
combined_tmr_remap_enrichment = dplyr::rename(combined_tmr_remap_enrichment, 
  "relative_enrichment" = enrichment, "query_name" = tf_name, "test_overlap_size" = estimate1)

# Rank results by enrichment
combined_tmr_remap_enrichment = mutate(group_by(combined_tmr_remap_enrichment, group), 
  rank = rank(-relative_enrichment), ranking_name = paste(query_name, rank(-relative_enrichment), sep = "_"))
combined_tmr_remap_enrichment = arrange(combined_tmr_remap_enrichment, desc(rank), group)
combined_tmr_remap_enrichment$ranking = factor(combined_tmr_remap_enrichment$ranking_name, unique(combined_tmr_remap_enrichment$ranking_name))

# Create plot showing enrichment of TF binding sites among TMRs
combined_tmr_remap_enrichment_plot = ggplot(combined_tmr_remap_enrichment, 
  aes(y = ranking, x = relative_enrichment, color = q_value, size = test_overlap_size)) + 
    geom_point() +
    geom_vline(xintercept = 1, linetype = "dashed") +
    scale_y_discrete(labels = function(x) gsub("_[0-9.]*$", "",  x), 
      expand = expansion(mult = c(0.05, 0.05))) +
    scale_colour_continuous(guide = guide_colorbar(order = 1, reverse = T)) + 
    labs(x = "Relative Enrichment", y = "Transcriptional Regulator", title = NULL, size = "Overlap Size") +
  facet_wrap("group", nrow = 2, dir = "v", scales = "free")
combined_tmr_remap_enrichment_plot = plotR::customize_ggplot_theme(combined_tmr_remap_enrichment_plot, 
  xlab = "Relative Enrichment", ylab = "Transcriptional Regulator", title = NULL, color_title = "Adjusted\np-value", 
  legend_title_size = 22, legend_text_size = 20, axis_title_size = 22) +
    theme(strip.background = element_blank(), axis.text.x = element_text(size = 14), 
      axis.text.y = element_text(size = 20), strip.text.x = element_text(size = 24), axis.title.y = element_text(size = 22))
combined_tmr_remap_enrichment_plot
ggsave(plot = combined_tmr_remap_enrichment_plot, "../figures/figure5C.pdf", width = 27, height = 15.2)