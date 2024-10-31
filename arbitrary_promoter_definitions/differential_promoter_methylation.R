# Calculate promoter methylation in CPGEA samples for different promoter definitions and perform differential methylation for different definitions

# Load required packages
library(dplyr)
library(cowplot)
library(GenomicRanges)
source("../auxillary_scripts/plotting_functions.R")
source("../auxillary_scripts/diff_meth_methylsig.R")
source("../example_differential_methylation_plot_functions.R")

### Test differential methylation using the different promoter definitions

# Get promoter definition list
promoter_definition_list = readRDS("promoter_definition_list.rds")

# Set names of promoters equal to transcript_id
promoter_definition_list = lapply(promoter_definition_list, function(x) setNames(x, x$transcript_id))

# Get path to CPGEA Meth RSE and convert to a methylation RSE
cpgea_meth_rse = HDF5Array::loadHDF5SummarizedExperiment("../auxillary_data/methylation_data/cpgea_wgbs_with_coverage_hg38/")

# Create a BPPARAM object
bpparam = BiocParallel::MulticoreParam(10)

# Perform differential analysis for different promoter definitions. Takes 25 minutes with 10 cores. 
system.time({promoter_diff_meth_results = lapply(promoter_definition_list, function(x) diff_meth_methylsig(meth_rse = cpgea_meth_rse, genomic_regions = x, 
  max_sites_per_chunk = floor(625000000/ncol(cpgea_meth_rse)), group_column = "condition", case = "Tumour", control = "Normal", BPPARAM = bpparam))})
saveRDS(promoter_diff_meth_results, "promoter_diff_meth_results.rds")

### Create a plot summarizing the differential methylation results for the different promoter definitions

# Load differential methylation results and convert into data.frames
promoter_diff_meth_results = readRDS("promoter_diff_meth_results.rds")
for(x in names(promoter_diff_meth_results)){
  promoter_diff_meth_results[[x]]$transcript_id = names(promoter_diff_meth_results[[x]])
}
promoter_diff_meth_results = lapply(promoter_diff_meth_results, data.frame)

# Combine promoter differential methylation results into a single table
promoter_diff_meth_results_df = bind_rows(promoter_diff_meth_results, .id = "definition")

# Convert meth_diff into a proportion
promoter_diff_meth_results_df$meth_diff = promoter_diff_meth_results_df$meth_diff/100

# Remove results with NA values
promoter_diff_meth_results_df = filter(promoter_diff_meth_results_df, !is.na(fdr))

# Denote whether promoters are hypermethylated, hypomethylated or unchanged 
promoter_diff_meth_results_df = mutate(promoter_diff_meth_results_df, 
  meth_change = case_when(
    fdr < 0.05 & meth_diff > 0 ~ "Hypermethylated",
    fdr < 0.05 & meth_diff < 0 ~ "Hypomethylated",
    fdr > 0.05 ~ "Unchanged"
    )
  )

# Convert meth_change to a factor
promoter_diff_meth_results_df$meth_change = factor(promoter_diff_meth_results_df$meth_change, levels = c("Hypomethylated", "Unchanged", "Hypermethylated"))

# Get proportion of hypermethylated, hypomethylated and unchanged promoters for each definition
promoter_diff_meth_results_df_summary = mutate(
  summarize(group_by(promoter_diff_meth_results_df, definition, meth_change), count = n()),
  freq = count/sum(count))

# Create a barplot of proportion of hypermethylated, hypomethylated and unchanged promoters for each definition
meth_change_proportions_barplot = ggplot(promoter_diff_meth_results_df_summary, 
  aes(y = count, x = definition, fill = meth_change, label = paste0(round(freq, 2)*100, "%"))) +
  geom_col(position = "dodge", color  = "black") +
  geom_text(mapping = aes(x = definition, y = count + 1000, group = meth_change), position = position_dodge(width = 0.9), size = 5)

# Adjust theme of barplot and save
meth_change_proportions_barplot = customize_ggplot_theme(meth_change_proportions_barplot, title = NULL, 
  xlab = "Promoter Definition", ylab = "Number of Promoters", fill_colors = c("#4B878BFF", "grey", "#D01C1FFF"), 
  scale_y = scale_y_continuous(limits = c(0, 60000), expand = expansion(mult = c(0, 0.05)), labels = scales::comma))
meth_change_proportions_barplot

### Create plots showing differential methylation results for example transcripts

# Create a data.frame with the coordinates of the different promoter definitions
promoter_definition_df = readRDS("promoter_definition_df.rds")

# Make promoter coordinates plot
promoter_region_plot = ggplot(promoter_definition_df, 
  aes(xmin = upstream, xmax = downstream, x = NULL, y = definition,  group = definition)) + 
  geom_linerange(linewidth = 12, position = position_dodge(0.06), color = RColorBrewer::brewer.pal(9, "YlGn")[c(2, 4, 6, 8, 9)]) +
  theme_classic() + 
  theme(plot.title = element_text(hjust = 0.5, size = 24), legend.text = element_text(size = 18),
    axis.title = element_text(size = 20), axis.text = element_text(size = 18))  +
  scale_x_continuous(limits = c(-5000, 5000), breaks = seq(-4000, 4000, 2000), expand = c(0.005, 0.005), labels = scales::comma) +
  geom_vline(xintercept = 0, linetype = "dashed") +
  labs(x = NULL, y = "Promoter\nDefinition")

# Add significance symbol to promoter_diff_meth_results
promoter_diff_meth_results$significance = ::sig_sym(promoter_diff_meth_results$fdr, symbol = "\u204E")

# Create a data.frame with the methylation values for all common_transcripts
promoter_definition_methylation_df = tidyr::pivot_wider(select(promoter_diff_meth_results, definition, transcript_id, meth_diff), names_from = definition, values_from = meth_diff)
promoter_definition_methylation_df = tibblecolumn_to_rownames(promoter_definition_methylation_df, "transcript_id")

# Make plots of CpG methylation change and promoter definition differential methylation for FLT1 and SLC5A8
flt1_cpg_meth_change_plot = plot_cpg_methylation_change(transcript = "ENST00000282397", 
  xlabel = expression("Distance to" ~ italic("FLT1") ~ "TSS (bp)"))
slc5a8_cpg_meth_change_plot = plot_cpg_methylation_change(transcript = "ENST00000536262", 
  xlabel = expression("Distance to" ~ italic("SLC5A8") ~ "TSS (bp)"))
flt1_promoters_plot = plot_promoter_methylation_change("ENST00000282397", 
  title = "*FLT1* Promoter<br>Methylation Change")
slc5a8_promoters_plot = plot_promoter_methylation_change("ENST00000536262", 
  title = "*SLC5A8* Promoter<br>Methylation Change")

# Combine plots into a single figure along with promoter definition plot and save
promoters_and_examples_plot_list = list(promoter_region_plot, NULL, 
  flt1_cpg_meth_change_plot, flt1_promoters_plot,
  slc5a8_cpg_meth_change_plot, slc5a8_promoters_plot)
promoters_and_examples_plot = plot_grid(plotlist = promoters_and_examples_plot_list, nrow = 3, ncol = 2, align = "hv", 
  rel_heights = c(1.5, 5.5, 5.5), rel_widths = c(3.5, 1), labels = c("A", "", "B", "", "C", ""))
promoters_and_examples_plot

# Add label "D" to barplot
meth_change_proportions_barplot = plot_grid(plotlist = list(meth_change_proportions_barplot), labels = "D")

# Add barplot to complete plot
figure2_plot_list = list(promoters_and_examples_plot, meth_change_proportions_barplot)
figure2_plot = plot_grid(plotlist = figure2_plot_list, nrow = 2, ncol = 1,
  rel_heights = c(0.65, 0.35))
ggsave(plot = figure2_plot, filename = "../figures/figure2.pdf", 
  width = 20.57, height = 32.72, device = cairo_pdf)

### Make hypermethylated/hypomethylated promoters overlaps plot 

# Get the hypermethylated and hypomethylated transcripts for each promoter definition
hypermethylated_promoters = lapply(promoter_diff_meth_results, function(x) 
  filter(x, fdr < 0.05, meth_diff > 0)$transcript_id)
hypomethylated_promoters = lapply(promoter_diff_meth_results, function(x) 
  filter(x, fdr < 0.05, meth_diff < 0)$transcript_id)

# Calculate the intersections for hypermethylated and hypomethylated transcripts for the different promoter definitions
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
  filename = "../figures//supp_figure5.pdf", height = 9, width = 9)
  
### Make upset diagrams for hypermethylated and hypomethylated promoters

hypermethylated_promoters = lapply(hypermethylated_promoters, function(x) intersect(x, cage_transcripts))
pdf(file = "../figures/hypermethylated_upset.pdf", width = 25, height = 14, bg = "white")
UpSetR::upset(UpSetR::fromList(hypermethylated_promoters), sets.bar.color = RColorBrewer::brewer.pal(9, "YlGn")[c(2, 4, 6, 8, 9)], sets = rev(names(promoter_diff_meth_results)), keep.order = T,
  decreasing = c(T, T), nsets = length(hypermethylated_promoters), nintersects = NA, text.scale = c(3.5, rep(3, 5)),  mainbar.y.label = "Intersection Size\n\n")
grid::grid.text("Overlaps of Significantly Hypermethylated Promoters\nfrom Different Promoter Definitions", x = 0.65, y=0.95, gp=grid::gpar(fontsize=30))
dev.off()

hypomethylated_promoters = lapply(hypomethylated_promoters, function(x) intersect(x, cage_transcripts))
pdf(file = "../figures/hypomethylated_upset.pdf", width = 25, height = 14, bg = "white")
UpSetR::upset(UpSetR::fromList(hypomethylated_promoters), sets.bar.color = RColorBrewer::brewer.pal(9, "YlGn")[c(2, 4, 6, 8, 9)], sets = rev(names(promoter_diff_meth_results)), keep.order = T,
  decreasing = c(T, T), nsets = length(hypomethylated_promoters), nintersects = NA, text.scale = c(3.5, rep(3, 5)),  mainbar.y.label = "Intersection Size\n\n", set_size.scale_max = 11000)
grid::grid.text("Overlaps of Significantly Hypomethylated Promoters\nfrom Different Promoter Definitions", x = 0.65, y=0.95, gp=grid::gpar(fontsize=30))
dev.off()