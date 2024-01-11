# Perform overrepresentation analysis of MSigDB Hallmark and KEGG pathways among TMR-associated genes

# Load required packages
library(dplyr)
library(enrichmentTests)
library(plotR)
library(ggpubr)

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

# Get genes associated with TMRs from each dataset
tmr_genes = lapply(tmr_list, function(x) unique(x$gene_name))

# Get all background protein-coding genes
background_genes = unique(readRDS("~/genomes/gencode/gencode_granges/pcg_transcript_tss_ranges_gencode_v38.rds")$gene_name)

# Load MSigDB gene sets
msigdb_gene_set_list = enrichmentTests::msigdb_gene_set_list

# Perform overrepresentation analysis for Hallmark pathways 
system.time({hallmark_enrichment_tmr_genes = lapply(tmr_genes, function(x) 
  fisher_test_apply(test = x, universe = background_genes, query_list = msigdb_gene_set_list$h, return_overlap = F))})
sapply(hallmark_enrichment_tmr_genes, function(x) sum(x$q_value < 0.05))

# Combine hallmark_enrichment_tmr_genes into a single table
combined_hallmark_results = bind_rows(c(hallmark_enrichment_tmr_genes), .id = "group")

# Filter for significant results
combined_hallmark_results = filter(combined_hallmark_results, q_value < 0.05)

# Rank the results separately for each dataset
combined_hallmark_results = mutate(group_by(combined_hallmark_results, group), 
  rank = rank(-relative_enrichment), ranking_name = paste0(query_name, rank(-relative_enrichment)))
combined_hallmark_results$group = factor(combined_hallmark_results$group, unique(combined_hallmark_results$group))
combined_hallmark_results = arrange(combined_hallmark_results, desc(rank), group)
combined_hallmark_results$ranking = factor(combined_hallmark_results$ranking_name, unique(combined_hallmark_results$ranking_name))

# Plot the Hallmark enrichment results
combined_hallmark_enrichment_plot = ggplot(combined_hallmark_results, 
  aes(y = ranking, x = relative_enrichment, color = q_value, size = test_overlap_size)) + 
    geom_point() + geom_vline(xintercept = 1, linetype = "dashed")
combined_hallmark_enrichment_plot = customize_ggplot_theme(combined_hallmark_enrichment_plot, xlab = "Relative Enrichment", ylab = "MSigDB Pathway", color_title = "Significance",
  facet = "group", facet_labels = c("Normal Prostate", "Prostate Tumours", "Prostate Metastases"), facet_scales = "free_y", strip_text_size = 20) + 
  theme(axis.text.y = element_text(size = 12), strip.background = element_blank()) +
  scale_y_discrete(labels = function(x) gsub("[0-9.]*$", "", gsub("_", " ", x))) +
  labs(size = "Overlap Size") +
  scale_colour_continuous(guide = guide_colorbar(order = 1, reverse = T))

# Perform overrepresentation analysis for KEGG pathways 
system.time({kegg_enrichment_tmr_genes = lapply(tmr_genes, function(x) 
  fisher_test_apply(test = x, universe = background_genes, query_list = msigdb_gene_set_list$`CP:KEGG`, return_overlap = F))})
sapply(kegg_enrichment_tmr_genes, function(x) sum(x$q_value < 0.05))

# Combine kegg_enrichment_tmr_genes into a single table
combined_kegg_results = bind_rows(c(kegg_enrichment_tmr_genes), .id = "group")

# Filter for significant results
combined_kegg_results = filter(combined_kegg_results, q_value < 0.05)
combined_kegg_results = mutate(group_by(combined_kegg_results, group), 
  rank = rank(-relative_enrichment), ranking_name = paste0(query_name, rank(-relative_enrichment)))

# Rank the results separately for each dataset
combined_kegg_results$group = factor(combined_kegg_results$group, unique(combined_kegg_results$group))
combined_kegg_results = arrange(combined_kegg_results, desc(rank), group)
combined_kegg_results$ranking = factor(combined_kegg_results$ranking_name, unique(combined_kegg_results$ranking_name))

combined_kegg_enrichment_plot = ggplot(combined_kegg_results, 
  aes(y = ranking, x = relative_enrichment, color = q_value, size = test_overlap_size)) + 
    geom_point() + geom_vline(xintercept = 1, linetype = "dashed")
combined_kegg_enrichment_plot = customize_ggplot_theme(combined_kegg_enrichment_plot, xlab = "Relative Enrichment", ylab = "KEGG Pathway", color_title = "Significance",
  facet = "group", facet_labels = c("Normal Prostate", "Prostate Tumours", "Prostate Metastases"), facet_scales = "free_y", strip_text_size = 20) + 
  theme(axis.text.y = element_text(size = 12), strip.background = element_blank()) +
  scale_y_discrete(labels = function(x) gsub("[0-9.]*$", "", gsub("_", " ", x))) +
  labs(size = "Overlap Size") +
  scale_colour_continuous(guide = guide_colorbar(order = 1, reverse = T))

# Combine Hallmark and KEGG enrichment results
combined_enrichment_plots = ggarrange(plotlist = list(combined_hallmark_enrichment_plot, combined_kegg_enrichment_plot), nrow = 2, align = "hv")
ggsave(plot = combined_enrichment_plots, "tmr_stats_plots/pathway_enrichment_plots.pdf", width = 27, height = 18)













###
tmr_diff_meth_genes = readRDS("tmr_methylation/tmr_diff_meth_genes.rds")

# Make a function which will take a list of genes and perform pathway enrichment
plot_tmr_gene_enrichment = function(gene_list, pathways, filter_groups = NULL, 
  facet_labels = c("Normal Prostate", "Prostate Tumours", "Prostate Metastases"), ylab = "Pathway", title = NULL){
  
  # Enrichment results
  enrichment_results = lapply(gene_list, function(x) 
    fisher_test_apply(test = x, universe = background_genes, query_list = pathways, return_overlap = F))

  # Combine enrichment_results into a single table
  combined_enrichment_results = bind_rows(c(enrichment_results), .id = "group")

  # Filter for significant results
  combined_enrichment_results = filter(combined_enrichment_results, q_value < 0.05)
  combined_enrichment_results = mutate(group_by(combined_enrichment_results, group), 
    rank = rank(-relative_enrichment), ranking_name = paste0(query_name, rank(-relative_enrichment)))

  # Rank the results separately for each dataset
  combined_enrichment_results$group = factor(combined_enrichment_results$group, unique(combined_enrichment_results$group))
  combined_enrichment_results = arrange(combined_enrichment_results, desc(rank), group)
  combined_enrichment_results$ranking = factor(combined_enrichment_results$ranking_name, unique(combined_enrichment_results$ranking_name))
  
  if(nrow(combined_enrichment_results) == 0){
    return(NULL)
  }
  
  if(!is.null(filter_groups)){
    combined_enrichment_results = dplyr::filter(combined_enrichment_results, group %in% filter_groups)
  }

  combined_enrichment_plot = ggplot(combined_enrichment_results, 
    aes(y = ranking, x = relative_enrichment, color = q_value, size = test_overlap_size)) + 
      geom_point() + geom_vline(xintercept = 1, linetype = "dashed")
  combined_enrichment_plot = customize_ggplot_theme(combined_enrichment_plot, title = title, xlab = "Relative Enrichment", ylab = ylab, color_title = "Adjusted\np-value",
    facet = "group", facet_labels = facet_labels, facet_scales = "free_y", strip_text_size = 20) + 
    theme(axis.text.y = element_text(size = 14), strip.background = element_blank()) +
    scale_y_discrete(labels = function(x) gsub("[0-9.]*$", "", gsub("_", " ", x))) +
    labs(size = "Overlap Size") +
    scale_colour_continuous(guide = guide_colorbar(order = 1, reverse = T))

  combined_enrichment_plot
  
}

system.time({kegg_enrichment_plots = lapply(tmr_diff_meth_genes, function(x) 
  plot_tmr_gene_enrichment(gene_list = x, pathways = msigdb_gene_set_list$`CP:KEGG`, ylab = "KEGG Pathway"))})

system.time({hallmark_enrichment_plots = lapply(tmr_diff_meth_genes, function(x) 
  plot_tmr_gene_enrichment(gene_list = x, pathways = msigdb_gene_set_list$h, ylab = "MSigDB Hallmark Pathway"))})

hypermethylated_negative_tmr_pathway_enrichment_plots = ggarrange(plotlist = list(kegg_enrichment_plots[[1]], hallmark_enrichment_plots[[1]]), nrow = 2, align = "hv", labels = c("A", "B"))
ggsave(plot = hypermethylated_negative_tmr_pathway_enrichment_plots, "tmr_stats_plots/hypermethylated_negative_tmr_pathway_enrichment_plots.pdf", width = 27, height = 18)

hypomethylated_negative_tmr_pathway_enrichment_plots = ggarrange(plotlist = list(kegg_enrichment_plots[[2]], hallmark_enrichment_plots[[2]]), nrow = 2, align = "hv", labels = c("A", "B"))
ggsave(plot = hypomethylated_negative_tmr_pathway_enrichment_plots, "tmr_stats_plots/hypomethylated_negative_tmr_pathway_enrichment_plots.pdf", width = 27, height = 18)

# Make presentation plots just for cancer
system.time({kegg_enrichment_plots_cancer = lapply(tmr_diff_meth_genes, function(x) 
  plot_tmr_gene_enrichment(gene_list = x, pathways = msigdb_gene_set_list$`CP:KEGG`, 
    filter_groups = c("cpgea_tumour", "mcrpc"), facet_labels = c("Prostate Tumours TMRs", "Prostate Metastases TMRs"), ylab = "KEGG Pathway"))})
system.time({msigdb_enrichment_plots_cancer = lapply(tmr_diff_meth_genes, function(x) 
  plot_tmr_gene_enrichment(gene_list = x, pathways = msigdb_gene_set_list$h, 
    filter_groups = c("cpgea_normal", "cpgea_tumour", "mcrpc"), facet_labels = c("Normal Prostate TMRs", "Prostate Tumours TMRs", "Prostate Metastases TMRs"), ylab = "MSigDB Hallmark Pathway"))})
ggsave(plot = kegg_enrichment_plots_cancer[[1]], "~/promoter_project/presentation_figures/hypermethylated_kegg_enrichment_plots_cancer.pdf", width = 27, height = 12)
ggsave(plot = msigdb_enrichment_plots_cancer[[2]], "~/promoter_project/presentation_figures/hypomethylated_msigdb_enrichment_plots_cancer.pdf", width = 27, height = 12)
supp_figure5_plots = ggarrange(plotlist = list(
  customize_ggplot_theme(kegg_enrichment_plots_cancer[[1]], title = "Hypermethylated Negative TMRs", xlab = "Relative Enrichment", ylab = "KEGG Pathway") + theme(strip.background = element_blank()), 
  customize_ggplot_theme(msigdb_enrichment_plots_cancer[[2]], title = "Hypomethylated Negative TMRs", xlab = "Relative Enrichment", ylab = "MSigDB Hallmark Pathway") + theme(strip.background = element_blank())),  
    nrow = 2, labels = c("A", "B"))
ggsave(plot = supp_figure5_plots, "tmr_stats_plots/supp_figure5_plots.pdf", width = 27, height = 18)
