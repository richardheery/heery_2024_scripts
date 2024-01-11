# Invesitgate methylation of TMRs in TCGA samples

# Load required packages
library(methodical)
library(foreach)
library(dplyr)
library(plotR)

# Load meth RSE for TCGA
tcga_meth_rse_hg19 = HDF5Array::loadHDF5SummarizedExperiment("~/tcga_methylation/tcga_meth_rse_hdf5_3_decimals_hg19/")

# Get hg19 to hg38 chain
hg19tohg38_chain = rtracklayer::import.chain("~/genomes/liftover_chain_files/hg19ToHg38.over.chain")

# Get CpG sites in hg38
hg38_cpgs = methodical:::cpg_genome_ranges_hg38

# Liftover tcga_meth_rse_hg19 to hg38
tcga_meth_rse_hg38 = liftover_meth_rse(meth_rse = tcga_meth_rse_hg19, chain = hg19tohg38_chain, 
  permitted_target_regions = hg38_cpgs)

# Get genome annotation for hg38
genome_annotation = methodical::genome_annotation_hg38

# Shorten names of exons and introns
genome_annotation$region_type[genome_annotation$region_type == "exon"] = "Exons"
genome_annotation$region_type[genome_annotation$region_type == "intron"] = "Introns"

# Remove " Region" from cluster_annotation$region_type
genome_annotation$region_type = gsub(" Region", "", genome_annotation$region_type)

# Remove lncRNA regions
genome_annotation = filter(genome_annotation, !grepl("lncRNA", region_type))

# Remove repeats, TSS, rRNA, small RNA, pseudogene and Protein-Coding from annotation
genome_annotation = plyranges::filter(genome_annotation, !region_type %in% 
    c("Alu", "CR1", "DNA Transposon", "L1", "L2", "LTR", "Low Complexity", "MIR", "SVA", 
      "Satellite", "Simple Repeat", "TSS", "rRNA", "Small RNA", "Pseudogene", "Protein-Coding"))

# Get repeat ranges for hg38
repeat_ranges = readRDS("~/genomes/repetitive_sequences/repeatmasker/repeatmasker_granges_ucsc.rds")

# Remove any features which overlap repeats
genome_annotation = subsetByOverlaps(genome_annotation, repeat_ranges, invert = T, ignore.strand = T)

# Get a list of TMRs
tmr_list = readRDS("../final_tmr_granges/final_tmr_5kb_list.rds")

# Create separate GRanges for prostate tumour negative and positive TMRs
negative_tmrs = unlist(GRangesList(tmr_list[c(3)]))
positive_tmrs = unlist(GRangesList(tmr_list[c(4)]))

# Replace metadata with one column giving the TMR direction
mcols(negative_tmrs) = data.frame(region_type = "Negative TMRs")
mcols(positive_tmrs) = data.frame(region_type = "Positive TMRs")

# Add TMRs to genome annotation
genome_annotation_with_tmrs = c(genome_annotation, negative_tmrs, positive_tmrs)
gc()

# Remove any regions not overlapping probes. Removes 85% of regions. At least 500 probes overlap each region type. 
genome_annotation_with_tmrs = subsetByOverlaps(genome_annotation_with_tmrs, rowRanges(tcga_meth_rse_hg38), ignore.strand = T)

# Create a list with each type of region
genome_annotation_with_tmrs_list = split(genome_annotation_with_tmrs, genome_annotation_with_tmrs$region_type)

# Get the mean methylation of probes for each genomic region. Took 50 minutes with 1 core locally
system.time({genomic_region_mean_probe_methylation = foreach(region = names(genome_annotation_with_tmrs_list), .packages = "methodical") %do% {
  
  message(paste("Starting", region))
  
  region_probes = subsetByOverlaps(tcga_meth_rse_hg38, genome_annotation_with_tmrs_list[[region]])
  system.time({region_methylation = colMeans(assay(region_probes), na.rm = T)})
  region_methylation
  
}})

# Add names to the results, combine them into a single table and save table
names(genomic_region_mean_probe_methylation) = names(genome_annotation_with_tmrs_list)
genomic_region_mean_probe_methylation = dplyr::bind_rows(genomic_region_mean_probe_methylation, .id = "region_type")
helpR::fwrite(genomic_region_mean_probe_methylation, "tmr_methylation_tables/tcga_genomic_region_mean_probe_methylation.tsv.gz")
genomic_region_mean_probe_methylation = data.frame(data.table::fread("tmr_methylation_tables/tcga_genomic_region_mean_probe_methylation.tsv.gz"), row.names = 1)

# Get the samples associated with each cancer type
tcga_samples_list = split(colnames(tcga_meth_rse_hg19), colData(tcga_meth_rse_hg19)$project)

# For each cancer, find the mean methylation change for each region
tcga_genomic_feature_mean_meth_change  = foreach(cancer = names(tcga_samples_list)) %do% {
  
  # Get samples for specified cancer type and get associated normal and tumour samples
  samples = tcga_samples_list[[cancer]]
  normal_samples = samples[endsWith(samples, "11")]
  tumour_samples = samples[endsWith(samples, "01")]
  
  # Get the mean methylation chnage for probes overlapping each region type and return
  region_mean_meth_change = rowMeans(select(genomic_region_mean_probe_methylation, all_of(tumour_samples)), na.rm = T) - 
    rowMeans(select(genomic_region_mean_probe_methylation, all_of(normal_samples)), na.rm = T)
  region_mean_meth_change
  
}

# Add names to results and combine into a single table
names(tcga_genomic_feature_mean_meth_change) = names(tcga_samples_list)
tcga_genomic_feature_mean_meth_change = bind_rows(tcga_genomic_feature_mean_meth_change, .id = "cancer")

# Remove cancers with missing values for all features
tcga_genomic_feature_mean_meth_change = tcga_genomic_feature_mean_meth_change[rowSums(is.na(tcga_genomic_feature_mean_meth_change)) < 10, ]

# Create a vector with the selected genomic features to plot
genomic_feature_plot_order = c("Negative TMRs", "Positive TMRs", "CpG Island", "Predicted Promoter", "Predicted Enhancer", 
  "Open Chromatin", "CTCF BS", "TF BS", "Exons", "Introns")

# Convert tcga_genomic_feature_mean_meth_change to long format
tcga_genomic_feature_mean_meth_change_long = tidyr::pivot_longer(tcga_genomic_feature_mean_meth_change, 
  cols = -cancer, names_to = "region_type", values_to = "mean_meth_change")

# Filter for selectd regions
tcga_genomic_feature_mean_meth_change_long_selected = filter(tcga_genomic_feature_mean_meth_change_long, 
  region_type %in% genomic_feature_plot_order)

# Convert region_type to a factor with levels in correct order
tcga_genomic_feature_mean_meth_change_long_selected$region_type = 
  factor(tcga_genomic_feature_mean_meth_change_long_selected$region_type, levels = genomic_feature_plot_order)

# Set colours for genomic features
feature_colors = c(colour_list$purple_and_gold_light, rev(RColorBrewer::brewer.pal(8, name = "BrBG")))

# Add a group for splitting the plot into two
tcga_genomic_feature_mean_meth_change_long_selected$group = 
  ceiling(as.numeric(factor(tcga_genomic_feature_mean_meth_change_long_selected$cancer))/
      (length(unique(tcga_genomic_feature_mean_meth_change_long_selected$cancer))/2))

# Reverse order of cancer type so they are in alphabetical order from top to bottom
tcga_genomic_feature_mean_meth_change_long_selected$cancer = 
  factor(tcga_genomic_feature_mean_meth_change_long_selected$cancer, levels = rev(unique(tcga_genomic_feature_mean_meth_change_long_selected$cancer)))

# Give TMRs a diamond shape
tcga_genomic_feature_mean_meth_change_long_selected$shape = 
  ifelse(grepl("TMR", tcga_genomic_feature_mean_meth_change_long_selected$region_type), 23, 21)

# Create a plot of the mean methylation change of probes for different genomic features in TCGA
tcga_feature_meth_change_plot = ggplot(tcga_genomic_feature_mean_meth_change_long_selected,
  aes(y = cancer, x = mean_meth_change, fill = region_type, color = region_type)) +
  geom_point(shape = tcga_genomic_feature_mean_meth_change_long_selected$shape,
   size = 6, position = position_jitter(seed = 2, height = 0.3), alpha = 0.7) +
  geom_vline(xintercept = 0, linetype = "dashed") +
  geom_hline(yintercept = seq(0.5, length(unique(tcga_genomic_feature_mean_meth_change_long_selected$cancer)), by = 1), 
    color = "gray", linewidth = .5, alpha = .5)
tcga_feature_meth_change_plot = 
  customize_ggplot_theme(tcga_feature_meth_change_plot, 
  ylab = "Cancer Type", xlab = "Mean Methylation Change", color_title = "Genomic Feature",
  show_legend = T, scale_x = scale_x_continuous(limits = c(-0.1, 0.12), expand = c(0, 0)), 
  fill_colors = feature_colors, colors = feature_colors, legend_text_size = 16, 
  facet = "group") + 
  theme(strip.background = element_blank(), strip.text = element_blank(),
    panel.grid.major.y = element_blank()) +
    guides(fill = "none")
tcga_feature_meth_change_plot
ggsave(plot = tcga_feature_meth_change_plot, filename = "tcga_feature_meth_change_plot.pdf", height = 9, width = 16)

# Create a vector with full cancer names
cancer_names = rep(c("Bladder", "Breast", "Cervical", "Cholangiocarcinoma", "Colon", "Esophageal",
  "Glioma", "Head and Neck", "Kidney Clear Cell", "Kidney Papillary Cell", "Liver",
  "Lung Adenocarcinoma", "Lung Squamous", "Pancreatic", "PCPG",
  "Prostate", "Rectal", "Sarcoma", "Skin Melanoma", "Stomach", "Thyroid", "Thymoma", "Endometrial"), each = 10)

tcga_genomic_feature_mean_meth_change_long_selected$full_cancer = factor(cancer_names, levels = rev(unique(cancer_names)))

tcga_feature_meth_change_plot = ggplot(tcga_genomic_feature_mean_meth_change_long_selected,
  aes(y = full_cancer, x = mean_meth_change, fill = region_type, color = region_type)) +
  geom_point(shape = tcga_genomic_feature_mean_meth_change_long_selected$shape,
   size = 6, position = position_jitter(seed = 2, height = 0.3), alpha = 0.7) +
  geom_vline(xintercept = 0, linetype = "dashed") +
  geom_hline(yintercept = seq(0.5, length(unique(tcga_genomic_feature_mean_meth_change_long_selected$cancer)), by = 1), 
    color = "gray", linewidth = .5, alpha = .5)
tcga_feature_meth_change_plot = 
  customize_ggplot_theme(tcga_feature_meth_change_plot, 
  ylab = "Cancer Type", xlab = "Mean Methylation Change", color_title = "TMRs and\nRegulatory Features",
  show_legend = T, scale_x = scale_x_continuous(limits = c(-0.1, 0.12), expand = c(0, 0)), 
  fill_colors = feature_colors, colors = feature_colors, legend_text_size = 16, 
  facet = "group") + 
  theme(strip.background = element_blank(), strip.text = element_blank(),
    panel.grid.major.y = element_blank()) +
    guides(fill = "none")
tcga_feature_meth_change_plot = ggpubr::ggarrange(tcga_feature_meth_change_plot, labels = "C")
tcga_feature_meth_change_plot
ggsave(plot = tcga_feature_meth_change_plot, filename = "~/promoter_project/presentation_figures/tcga_feature_meth_change_plot.1.pdf", height = 9, width = 16)

### Evalaute region methylation in TCGA WGBS data

# Load meth RSE for TCGA WGBS data
tcga_wgbs_meth_rse = HDF5Array::loadHDF5SummarizedExperiment("~/wgbs/tcga_wgbs/tcga_wgbs_hg38/")

# Get the mean methylation of CpGs for each genomic region. Took 10 minutes with 1 core locally
system.time({tcga_wgbs_genomic_region_mean_meth = foreach(region = names(genome_annotation_with_tmrs_list), .packages = "methodical") %do% {
  
  message(paste("Starting", region))
  
  region_cpgs = subsetByOverlaps(tcga_wgbs_meth_rse, genome_annotation_with_tmrs_list[[region]])
  system.time({region_methylation = colMeans(assay(region_cpgs), na.rm = T)})
  region_methylation
  
}})

# Add names to results and combine into a single table
names(tcga_wgbs_genomic_region_mean_meth) = names(genome_annotation_with_tmrs_list)
tcga_wgbs_genomic_region_mean_meth = dplyr::bind_rows(tcga_wgbs_genomic_region_mean_meth, .id = "region_type")
helpR::fwrite(tcga_wgbs_genomic_region_mean_meth, "tmr_methylation_tables/tcga_wgbs_genomic_region_mean_meth.tsv.gz")
tcga_wgbs_genomic_region_mean_meth = data.frame(data.table::fread("tmr_methylation_tables/tcga_wgbs_genomic_region_mean_meth.tsv.gz"), row.names = 1)

# Identify submitters with matching tumour and normal samples. There is one submitter each for BLCA, BRCA, COAD, LUAD, LUSC, READ, STAD and UCEC
normal_samples = rownames(colData(tcga_wgbs_meth_rse))[which(colData(tcga_wgbs_meth_rse)$sample_type == 11)]
matching_tumour_samples = gsub("_11", "_1", normal_samples)

# Find mean methylation change for all regions 
tcga_wgbs_genomic_region_mean_meth_tumour = select(tcga_wgbs_genomic_region_mean_meth, all_of(matching_tumour_samples))
tcga_wgbs_genomic_region_mean_meth_normal = select(tcga_wgbs_genomic_region_mean_meth, all_of(normal_samples))
tcga_wgbs_genomic_region_mean_meth_change = tcga_wgbs_genomic_region_mean_meth_tumour - tcga_wgbs_genomic_region_mean_meth_normal

# Change row names to a column and convert to long format
tcga_wgbs_genomic_region_mean_meth_change = tibble::rownames_to_column(tcga_wgbs_genomic_region_mean_meth_change, "region_type")
tcga_wgbs_genomic_region_mean_meth_change = tidyr::pivot_longer(tcga_wgbs_genomic_region_mean_meth_change, 
  cols = -region_type, names_to = "sample_name", values_to = "mean_meth_change")

# Add tumour type
tcga_wgbs_genomic_region_mean_meth_change$cancer = colData(tcga_wgbs_meth_rse)[gsub("_1", "_01", tcga_wgbs_genomic_region_mean_meth_change$sample_name), ]$project

# Give TMRs a diamond shape
tcga_wgbs_genomic_region_mean_meth_change$shape = 
  ifelse(grepl("TMR", tcga_wgbs_genomic_region_mean_meth_change$region_type), 23, 21)

# Convert region_type to a factor with levels in correct order
tcga_wgbs_genomic_region_mean_meth_change$region_type = 
  factor(tcga_wgbs_genomic_region_mean_meth_change$region_type, levels = genomic_feature_plot_order)

# Reverse order of cancer type so they are in alphabetical order from top to bottom
tcga_wgbs_genomic_region_mean_meth_change$cancer = 
  factor(tcga_wgbs_genomic_region_mean_meth_change$cancer, levels = rev(unique(tcga_wgbs_genomic_region_mean_meth_change$cancer)))
levels(tcga_wgbs_genomic_region_mean_meth_change$cancer) = 
  c("Endometrial", "Stomach", "Rectal", "Lung Squamous", "Lung Adenocarcinoma", "Colon", "Breast", "Bladder")

# Create a plot of the mean methylation change of probes for different genomic features in TCGA
tcga_wgbs_feature_meth_change_plot = ggplot(tcga_wgbs_genomic_region_mean_meth_change,
  aes(y = cancer, x = mean_meth_change, fill = region_type, color = region_type)) +
  geom_point(shape = tcga_wgbs_genomic_region_mean_meth_change$shape,
   size = 6, position = position_jitter(seed = 2, height = 0.3), alpha = 0.7) +
  geom_vline(xintercept = 0, linetype = "dashed") +
  geom_hline(yintercept = seq(0.5, length(unique(tcga_genomic_feature_mean_meth_change_long_selected$cancer)), by = 1), 
    color = "gray", linewidth = .5, alpha = .5)
tcga_wgbs_feature_meth_change_plot = 
  customize_ggplot_theme(tcga_wgbs_feature_meth_change_plot, 
  ylab = "Cancer Type", xlab = "Mean Methylation Change", color_title = "Genomic Feature",
  show_legend = T, scale_x = scale_x_continuous(limits = c(-0.15, 0.3), expand = c(0, 0)), 
  fill_colors = feature_colors, colors = feature_colors, legend_text_size = 16) + 
  theme(strip.background = element_blank(), strip.text = element_blank(),
    panel.grid.major.y = element_blank()) +
    guides(fill = "none")
tcga_wgbs_feature_meth_change_plot
ggsave(plot = tcga_wgbs_feature_meth_change_plot, filename = "tcga_wgbs_feature_meth_change_plot.pdf", height = 9, width = 16)
