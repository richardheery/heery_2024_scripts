# Calculate promoter methylation in CPGEA samples for different promoter definitions and perform differential methylation for different definitions

# Load required packages
library(dplyr)
library(methodical)
library(ggplot2)
library(doParallel)
library(plotR)
source("~/r_functions/table_tests.R")

# Get promoter definition list
promoter_definition_list = readRDS("promoter_definition_list.rds")

# Get names of MANE transcripts
mane_transcripts = gsub("\\.[0-9]*", "", readRDS("~/genomes/gencode/gencode_granges/gencode_v38_mane_transcript_ids.rds"))

# Get path to CPGEA methrix and convert a methylation RSE
cpgea_rse = HDF5Array::loadHDF5SummarizedExperiment("~/wgbs/cpgea/wgbs/cpgea_meth_rse_hdf5")

# Get methylation values for all promoter methylation definitions and save to promoter_definition_methylation_tables. Took about 1 hour. 
system.time({for(definition in names(promoter_definition_list)){
  promoter_definition_methylation = methodicalFinal::summarize_region_methylation(
    meth_rse = cpgea_rse, genomic_regions_chunk_size = 5000, 
    genomic_regions = promoter_definition_list[[definition]], 
    genomic_regions_names = promoter_definition_list[[definition]]$transcript_id
    )
  fwrite(promoter_definition_methylation, 
    filename = sprintf("promoter_definition_methylation_tables/%s_definition_methylation_table.tsv.gz", definition), 
    sep = "\t", row.names = F, na = "NA", quote = F)
}})

# Get paths to all promoter methylation definition tables
promoter_definition_methylation_tables = list.files("promoter_definition_methylation_tables", full.names = T, pattern = "_methylation_table.tsv.gz")
names(promoter_definition_methylation_tables) = gsub("_definition_methylation_table.tsv.gz", "", basename(promoter_definition_methylation_tables))

# Read in tables as a list and save
promoter_definition_methylation_tables = lapply(promoter_definition_methylation_tables, function(x)
  data.frame(data.table::fread(x), row.names = 1))
saveRDS(promoter_definition_methylation_tables, "promoter_definition_methylation_tables/promoter_definition_methylation_tables.rds")
promoter_definition_methylation_tables = readRDS("promoter_definition_methylation_tables/promoter_definition_methylation_tables.rds")

# For each promoter definition perform differential methylation testing using paired t-tests
system.time({for(definition in names(promoter_definition_methylation_tables)){
  print(paste("starting", definition))
  methylation_table = promoter_definition_methylation_tables[[definition]]
  differential_methylation_results =
    table_tests_unmatching(test_table1 = methylation_table, group1_pattern = "T", group2_pattern = "N", paired = T, gene_names = row.names(methylation_table))
  output = paste0(definition, "_definition_differential_methylation_results.tsv.gz")
  data.table::fwrite(differential_methylation_results, paste0("promoter_differential_methylation_results/", output), sep = "\t", row.names = F, quote = F)
}})

# Get paths to differential methylation results
promoter_differential_methylation_tables = list.files("~/promoter_project/arbitrary_promoter_definitions/promoter_differential_methylation_results", 
  pattern = "tsv.gz", full.names = T)

# Add names to tables and convert to letters and put in order
names(promoter_differential_methylation_tables) = 
  gsub("_definition_differential_methylation_results.tsv.gz", "", basename(promoter_differential_methylation_tables))

# Create a list with differential methylation results for all promoters
promoter_diff_methylation_results = lapply(promoter_differential_methylation_tables, function(x)
  data.table::fread(x))
saveRDS(promoter_diff_methylation_results, "promoter_differential_methylation_results/promoter_diff_methylation_results.rds")
promoter_diff_methylation_results = readRDS("promoter_differential_methylation_results/promoter_diff_methylation_results.rds")

# Make a list with the methylation change for each promoter using each definition
promoter_definition_methylation_change = lapply(promoter_diff_methylation_results, function(x) 
  setNames(x$mean_change, x$feature))
saveRDS(promoter_definition_methylation_change, "promoter_definition_methylation_change.rds")

# Get the hypermethylated and hypomethylated transcripts for each promoter definition
hypermethylated_promoters = lapply(promoter_diff_methylation_results, function(x) 
  filter(x, q_value < 0.05, mean_change > 0)$gene)
hypomethylated_promoters = lapply(promoter_diff_methylation_results, function(x) 
  filter(x, q_value < 0.05, mean_change < 0)$gene)
saveRDS(list(hypermethylated_promoters = hypermethylated_promoters, hypomethylated_promoters = hypomethylated_promoters),
  "promoter_differential_methylation_results/differentially_methylated_promoters.rds")

# Calculate the intersections for hypermethylated and hypomethylated transcripts for the different promoter defintions
hypermethylated_hypomethylated_intersections = lapply(hypermethylated_promoters, function(x) 
  sapply(hypomethylated_promoters, function(y) length(intersect(x, y))))

# Convert result into a data.frame where columns are hypermethylated and rows are hypomethylated promoters
hypermethylated_hypomethylated_intersections = data.frame(hypermethylated_hypomethylated_intersections)

# Set the diagonal to NA
diag(hypermethylated_hypomethylated_intersections) = NA

# Convert rownames to a column called hypomethylated
hypermethylated_hypomethylated_intersections = tibble::rownames_to_column(hypermethylated_hypomethylated_intersections, "hypomethylated")

# Convert the results into long format
hypermethylated_hypomethylated_intersections = tidyr::pivot_longer(hypermethylated_hypomethylated_intersections, cols = -hypomethylated, names_to = "hypermethylated")

# Reverse the factor levels for hypermethylated
hypermethylated_hypomethylated_intersections$hypermethylated = forcats::fct_rev(hypermethylated_hypomethylated_intersections$hypermethylated)

# Create a plot with the size of the intersections
hypermethylated_hypomethylated_intersections_plot = 
  ggplot(hypermethylated_hypomethylated_intersections, 
    aes(x = hypomethylated, y = hypermethylated, size = value, fill = value, label = scales::comma(value))) +
  geom_point(shape = 21, color = "black") +
  geom_text(color = "white", size = 6) +
  theme_classic() +
  scale_radius(range = c(20, 40)) +
  labs(x = "Hypomethylated Promoters", y = "Hypermethylated Promoters", title = "Overlap of Differentially Methylated\nPromoters using Different Definitions") +
  theme(plot.title = element_text(hjust = 0.5, size = 20), axis.title = element_text(size = 20), 
    axis.text = element_text(size = 14), legend.position = "None")
ggsave(plot = hypermethylated_hypomethylated_intersections_plot, 
  filename = "promoter_definition_plots/differential_methylation_plots/hypermethylated_hypomethylated_intersections.pdf", height = 9, width = 9)
  
# Make upset diagram of hypermethylated promoters
pdf(file = "promoter_definition_plots/differential_methylation_plots/hypermethylated_upset.pdf", width = 25, height = 14, bg = "white")
UpSetR::upset(UpSetR::fromList(hypermethylated_promoters), sets.bar.color = RColorBrewer::brewer.pal(9, "YlGn")[c(2, 4, 6, 8, 9)], sets = rev(names(promoter_diff_methylation_results)), keep.order = T,
  decreasing = c(T, T), nsets = length(hypermethylated_promoters), nintersects = NA, text.scale = c(3.5, rep(3, 5)),  mainbar.y.label = "Intersection Size\n\n")
grid::grid.text("Overlaps of Significantly Hypermethylated Promoters\nfrom Different Promoter Definitions", x = 0.65, y=0.95, gp=grid::gpar(fontsize=30))
dev.off()

# Make upset diagram of hypomethylated promoters
pdf(file = "promoter_definition_plots/differential_methylation_plots/hypomethylated_upset.pdf", width = 25, height = 14, bg = "white")
UpSetR::upset(UpSetR::fromList(hypomethylated_promoters), sets.bar.color = RColorBrewer::brewer.pal(9, "YlGn")[c(2, 4, 6, 8, 9)], sets = rev(names(promoter_diff_methylation_results)), keep.order = T,
  decreasing = c(T, T), nsets = length(hypomethylated_promoters), nintersects = NA, text.scale = c(3.5, rep(3, 5)),  mainbar.y.label = "Intersection Size\n\n", set_size.scale_max = 42000)
grid::grid.text("Overlaps of Significantly Hypomethylated Promoters\nfrom Different Promoter Definitions", x = 0.65, y=0.95, gp=grid::gpar(fontsize=30))
dev.off()

### Make barplot of number of hypermethylated, hypomethylated and unchanged promoters for each definition

# Get differential methylation results for the different promoter definitions
promoter_diff_methylation_results = readRDS("promoter_differential_methylation_results/promoter_diff_methylation_results.rds")

# Combine the differential promoter methylation results into a single table
all_diff_meth_results = bind_rows(promoter_diff_methylation_results, .id = "group")

# Remove results with NA values
all_diff_meth_results = filter(all_diff_meth_results, !is.na(q_value))

# Denote whether promoters are hypermethylated, hypomethylated or unchanged 
all_diff_meth_results = mutate(all_diff_meth_results, 
  meth_change = case_when(
    q_value < 0.05 & mean_change > 0 ~ "Hypermethylated",
    q_value < 0.05 & mean_change < 0 ~ "Hypomethylated",
    q_value > 0.05 ~ "Unchanged"
    )
  )

# Convert meth_change to a factor
all_diff_meth_results$meth_change = factor(all_diff_meth_results$meth_change, levels = c("Hypomethylated", "Unchanged", "Hypermethylated"))

# Get proportion of hypermethylated, hypomethylated and unchanged promoters for each definition
all_diff_meth_results_summary = mutate(
  summarize(group_by(all_diff_meth_results, group, meth_change), count = n()),
  freq = count/sum(count))

# Create a barplot of proportion of hypermethylated, hypomethylated and unchanged promoters for each definition
meth_change_proportions_barplot = ggplot(all_diff_meth_results_summary, 
  aes(y = count, x = group, fill = meth_change, label = paste0(round(freq, 2)*100, "%"))) +
  geom_col(position = "dodge", color  = "black") +
  geom_text(mapping = aes(x = group, y = count + 1000, group = meth_change), position = position_dodge(width = 0.9), size = 5)

# Adjust theme of barplot and save
meth_change_proportions_barplot = customize_ggplot_theme(meth_change_proportions_barplot, title = NULL, 
  xlab = "Promoter Definition", ylab = "Number of Promoters", fill_colors = c("#4B878BFF", "grey", "#D01C1FFF"), 
  scale_y = scale_y_continuous(limits = c(0, 60000), expand = expansion(mult = c(0, 0.05)), labels = scales::comma))
meth_change_proportions_barplot
ggsave(plot = meth_change_proportions_barplot, width = 16, height = 9,
  "promoter_definition_plots/differential_methylation_plots/meth_change_proportions_barplot.pdf")
ggsave(plot = meth_change_proportions_barplot + labs(title = "Differentially Methylated Promoters in Prostate Cancer"), width = 16, height = 9,
  "~/promoter_project/presentation_figures/meth_change_proportions_barplot.pdf")
saveRDS(meth_change_proportions_barplot, "promoter_definition_plots/differential_methylation_plots/meth_change_proportions_barplot.rds")

meth_change_proportions_barplot_pres = ggplot(filter(all_diff_meth_results_summary, group %in% c("A", "E")),
  aes(y = count, x = group, fill = meth_change, label = paste0(round(freq, 2)*100, "%"))) +
  geom_col(position = "dodge", color  = "black") +
  geom_text(mapping = aes(x = group, y = count + 1000, group = meth_change), position = position_dodge(width = 0.9), size = 5)
# Adjust theme of barplot and save
meth_change_proportions_barplot_pres = customize_ggplot_theme(meth_change_proportions_barplot_pres, title = NULL, 
  xlab = "Promoter Definition", ylab = "Number of Promoters", fill_colors = c("#4B878BFF", "grey", "#D01C1FFF"), 
  fill_labels = c("Decreased Methylation", "Unchanged", "Increased Methylation"),
  scale_y = scale_y_continuous(limits = c(0, 60000), expand = expansion(mult = c(0, 0.05)), labels = scales::comma))
ggsave(plot = meth_change_proportions_barplot_pres + labs(title = "Differentially Methylated Promoters in Prostate Cancer"), width = 16, height = 9,
  "~/promoter_project/presentation_figures/meth_change_proportions_barplot_a_e.pdf")

### Make plots just with MANE transcripts

# Get the hypermethylated and hypomethylated transcripts for each promoter definition
hypermethylated_promoters_mane = lapply(promoter_diff_methylation_results, function(x) 
  intersect(filter(x, q_value < 0.05, mean_change > 0)$gene, mane_transcripts))
hypomethylated_promoters_mane = lapply(promoter_diff_methylation_results, function(x) 
  intersect(filter(x, q_value < 0.05, mean_change < 0)$gene, mane_transcripts))

# Calculate the intersections for hypermethylated and hypomethylated transcripts for the different promoter defintions
hypermethylated_hypomethylated_intersections_mane = lapply(hypermethylated_promoters_mane, function(x) 
  sapply(hypomethylated_promoters_mane, function(y) length(intersect(x, y))))

# Convert result into a data.frame where columns are hypermethylated and rows are hypomethylated promoters
hypermethylated_hypomethylated_intersections_mane = data.frame(hypermethylated_hypomethylated_intersections_mane)

# Set the diagonal to NA
diag(hypermethylated_hypomethylated_intersections_mane) = NA

# Convert rownames to a column called hypomethylated
hypermethylated_hypomethylated_intersections_mane = tibble::rownames_to_column(hypermethylated_hypomethylated_intersections_mane, "hypomethylated")

# Convert the results into long format
hypermethylated_hypomethylated_intersections_mane = tidyr::pivot_longer(hypermethylated_hypomethylated_intersections_mane, cols = -hypomethylated, names_to = "hypermethylated")

# Reverse the factor levels for hypermethylated
hypermethylated_hypomethylated_intersections_mane$hypermethylated = forcats::fct_rev(hypermethylated_hypomethylated_intersections_mane$hypermethylated)

# Create a plot with the size of the intersections
hypermethylated_hypomethylated_intersections_mane_plot = 
  ggplot(hypermethylated_hypomethylated_intersections_mane, 
    aes(x = hypomethylated, y = hypermethylated, size = value, fill = value, label = scales::comma(value))) +
  geom_point(shape = 21, color = "black") +
  geom_text(color = "white", size = 6) +
  theme_classic() +
  scale_radius(range = c(20, 40)) +
  labs(x = "Hypomethylated Promoters", y = "Hypermethylated Promoters", title = "Overlap of Differentially Methylated\nPromoters using Different Definitions") +
  theme(plot.title = element_text(hjust = 0.5, size = 20), axis.title = element_text(size = 20), 
    axis.text = element_text(size = 14), legend.position = "None")
ggsave(plot = hypermethylated_hypomethylated_intersections_mane_plot, 
  filename = "promoter_definition_plots/mane_plots/hypermethylated_hypomethylated_intersections_mane.pdf", height = 9, width = 9)

# Make upset diagram of hypermethylated promoters
pdf(file = "promoter_definition_plots/mane_plots/hypermethylated_upset_mane.pdf", width = 25, height = 14, bg = "white")
UpSetR::upset(UpSetR::fromList(hypermethylated_promoters_mane), sets.bar.color = RColorBrewer::brewer.pal(9, "YlGn")[c(2, 4, 6, 8, 9)], sets = rev(names(promoter_diff_methylation_results)), keep.order = T,
  decreasing = c(T, T), nsets = length(hypermethylated_promoters_mane), nintersects = NA, text.scale = c(3.5, rep(3, 5)),  mainbar.y.label = "Intersection Size\n\n")
grid::grid.text("Overlaps of Significantly Hypermethylated Promoters\nfrom Different Promoter Definitions", x = 0.65, y=0.95, gp=grid::gpar(fontsize=30))
dev.off()

# Make upset diagram of hypomethylated promoters
pdf(file = "promoter_definition_plots/mane_plots/hypomethylated_upset_mane.pdf", width = 25, height = 14, bg = "white")
UpSetR::upset(UpSetR::fromList(hypomethylated_promoters_mane), sets.bar.color = RColorBrewer::brewer.pal(9, "YlGn")[c(2, 4, 6, 8, 9)], sets = rev(names(promoter_diff_methylation_results)), keep.order = T,
  decreasing = c(T, T), nsets = length(hypomethylated_promoters_mane), nintersects = NA, text.scale = c(3.5, rep(3, 5)),  mainbar.y.label = "Intersection Size\n\n")
grid::grid.text("Overlaps of Significantly hypomethylated Promoters\nfrom Different Promoter Definitions", x = 0.65, y=0.95, gp=grid::gpar(fontsize=30))
dev.off()

# Convert meth_change to a factor
all_diff_meth_results_mane = filter(all_diff_meth_results, feature %in% mane_transcripts)

# Get proportion of hypermethylated, hypomethylated and unchanged promoters for each definition
all_diff_meth_results_mane_summary = mutate(
  summarize(group_by(all_diff_meth_results_mane, group, meth_change), count = n()),
  freq = count/sum(count))

# Create a barplot of proportion of hypermethylated, hypomethylated and unchanged promoters for each definition
meth_change_proportions_barplot_mane = ggplot(all_diff_meth_results_mane_summary, 
  aes(y = count, x = group, fill = meth_change, label = paste0(round(freq, 2)*100, "%"))) +
  geom_col(position = "dodge", color  = "black") +
  geom_text(mapping = aes(x = group, y = count + 1000, group = meth_change), position = position_dodge(width = 0.9), size = 5)

# Adjust theme of barplot and save
meth_change_proportions_barplot_mane = customize_ggplot_theme(meth_change_proportions_barplot_mane, title = NULL, 
  xlab = "Promoter Definition", ylab = "Number of Promoters", fill_colors = c("#4B878BFF", "grey", "#D01C1FFF"), 
  scale_y = scale_y_continuous(limits = c(0, 15000), expand = expansion(mult = c(0, 0.05)), labels = scales::comma))
meth_change_proportions_barplot_mane
ggsave(plot = meth_change_proportions_barplot_mane + labs(title = "Differentially Methylated Promoters in Prostate Cancer"), width = 16, height = 9,
  "promoter_definition_plots/mane_plots/meth_change_proportions_barplot_mane.pdf")
