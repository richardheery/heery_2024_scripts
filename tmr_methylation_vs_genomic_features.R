# Compare methylation change at TMRs to other genomic features. 

# Load required packages
library(methodical)
library(dplyr)
library(plotR)

# Get methylation values from CPGEA 
cpgea_rse = HDF5Array::loadHDF5SummarizedExperiment("~/mounts/local_mount/wgbs/cpgea/wgbs/cpgea_meth_rse_hdf5/")

# Get genome annotation for hg38
genome_annotation = methodicalFinal::genome_annotation_hg38

# Shorten names of exons and introns
genome_annotation$region_type[genome_annotation$region_type == "exon"] = "Exons"
genome_annotation$region_type[genome_annotation$region_type == "intron"] = "Introns"

# Remove " Region" from cluster_annotation$region_type
genome_annotation$region_type = gsub(" Region", "", genome_annotation$region_type)

# Remove lncRNA regions
genome_annotation = filter(genome_annotation, !grepl("lncRNA", region_type))

# Remove repeats and TSS from annotation
genome_annotation = plyranges::filter(genome_annotation, !region_type %in% 
    c("Alu", "CR1", "DNA Transposon", "L1", "L2", "LTR", "Low Complexity", "MIR", "SVA", "Satellite", "Simple Repeat", "TSS", "rRNA"))

# Get repeat ranges for hg38
repeat_ranges = readRDS("~/mounts/local_mount/genomes/repetitive_sequences/repeatmasker/repeatmasker_granges_ucsc.rds")

# Remove any features which overlap repeats
genome_annotation = subsetByOverlaps(genome_annotation, repeat_ranges, invert = T, ignore.strand = T)

# Get a list of TMRs
tmr_list = readRDS("../final_tmr_granges/final_tmr_5kb_list.rds")

# Rename TMRs
names(tmr_list) = c("Normal Prostate TMRs -", "Normal Prostate TMRs +", "Prostate Tumour TMRs -", 
  "Prostate Tumour TMRs +", "Prostate Metastasis TMRs -", "Prostate Metastasis TMRs +")

# Combine into a single GRanges
tmrs_gr = unlist(GRangesList(tmr_list))

# Remove metadata columns
mcols(tmrs_gr) = NULL

# Create a column with TMR group
tmrs_gr$region_type = names(tmrs_gr)
names(tmrs_gr) = NULL

# Add TMRs to genome annotation
genome_annotation_with_tmrs = c(genome_annotation, tmrs_gr)
gc()

# Find the number of probes overlapping each region
sort(sapply(split(genome_annotation_with_tmrs, genome_annotation_with_tmrs$region_type), function(x)
  length(subsetByOverlaps(rowRanges(tcga_meth_rse_hg38), x))))

# Get methylation values for genomic features and TMRs. Took 3 minutes with 10 cores. 
system.time({genomic_features_with_tmrs_methylation = 
  summarize_region_methylation(meth_rse = cpgea_rse, genomic_regions = genome_annotation_with_tmrs, keep_metadata_cols = T, n_chunks_parallel = 10)})

# Save table
helpR::fwrite(genomic_features_with_tmrs_methylation, "tmr_methylation_tables/genomic_features_with_tmrs_methylation.tsv.gz")

# Load methylation for genomic features and TMRs
genomic_features_with_tmrs_methylation = data.table::fread("tmr_methylation_tables/genomic_features_with_tmrs_methylation.tsv.gz")

# Drop region_name 
genomic_features_with_tmrs_methylation$region_name = NULL
gc()

# Split by region_type
genomic_features_with_tmrs_methylation = split(
  select(genomic_features_with_tmrs_methylation, -region_type), genomic_features_with_tmrs_methylation$region_type)
gc()

# Get the mean methylation change for each region
genomic_feature_mean_methylation_change = lapply(genomic_features_with_tmrs_methylation, function(x)
  rowMeans(select(x, starts_with("T")) - select(x, starts_with("N")), na.rm = T))

# Convert genomic_feature_mean_methylation_change into a single table and save
genomic_feature_mean_methylation_change = data.frame(
  region_type = rep(names(genomic_feature_mean_methylation_change), times = lengths(genomic_feature_mean_methylation_change)),
  values = unlist(genomic_feature_mean_methylation_change), row.names = NULL)
helpR::fwrite(genomic_feature_mean_methylation_change, "tmr_methylation_tables/genomic_feature_mean_methylation_change.tsv.gz")

# Load genomic_feature_mean_methylation_change
genomic_feature_mean_methylation_change = data.table::fread("tmr_methylation_tables/genomic_feature_mean_methylation_change.tsv.gz")

# Create a vector with the selected genomic features to plot
selected_genomic_features = c("CpG Island", "Predicted Promoter", "Predicted Enhancer", "Open Chromatin", 
  "CTCF BS", "TF BS", "Exons", "Introns")

# Filter genomic_feature_mean_methylation_change for selected features
genomic_feature_mean_methylation_change_selected = filter(genomic_feature_mean_methylation_change, 
  region_type %in% c(selected_genomic_features, c("Prostate Tumour TMRs -", "Prostate Tumour TMRs +")))

# Put region_type in desired order
genomic_feature_mean_methylation_change_selected$region_type = 
  factor(genomic_feature_mean_methylation_change_selected$region_type, 
    levels = rev(c(names(tmr_list)[c(1, 3, 5, 2, 4, 6)], selected_genomic_features)))

# Set colours for genomic features
feature_colors = rev(c(rep(colour_list$purple_and_gold_light, each = 1), 
  rev(RColorBrewer::brewer.pal(8, name = "BrBG"))))

# Randomly select 250 regions for each region type for plotting points
set.seed(123)
sample_regions = slice_sample(group_by(genomic_feature_mean_methylation_change_selected, region_type), n = 250)

# Create boxplots of mean methylation change for genomic features and TMRs
methylation_change_boxplots = ggplot(genomic_feature_mean_methylation_change_selected, 
  aes(x = region_type, y = values, fill = region_type)) +
  geom_boxplot(outlier.shape = NA) +
  geom_jitter(data = sample_regions, alpha = 1, aes(fill = region_type), shape = 21) +
  geom_hline(yintercept = 0, linetype = "dashed") +
  coord_flip()
methylation_change_boxplots = customize_ggplot_theme(methylation_change_boxplots, 
  ylab = "Mean Methylation Change", xlab = "TMRs and Regulatory Features", show_legend = F, 
  scale_y = scale_y_continuous(limits = c(-0.5, 0.5)), 
  fill_colors = feature_colors)
methylation_change_boxplots = ggpubr::ggarrange(methylation_change_boxplots, labels = "B")
methylation_change_boxplots
ggsave(plot = methylation_change_boxplots, "methylation_change_boxplots.1.pdf", width = 8, height = 9)

# 
genomic_feature_mean_methylation_change_selected_pres = filter(genomic_feature_mean_methylation_change_selected, 
  !grepl("Normal|Metastasis", region_type))
sample_regions_pres = filter(sample_regions, !grepl("Normal|Metastasis", region_type))
methylation_change_boxplots = ggplot(genomic_feature_mean_methylation_change_selected_pres, 
  aes(x = region_type, y = values, fill = region_type)) +
  geom_boxplot(outlier.shape = NA) +
  geom_jitter(data = sample_regions_pres, alpha = 1, aes(fill = region_type), shape = 21) +
  geom_hline(yintercept = 0, linetype = "dashed") +
  coord_flip()
methylation_change_boxplots = customize_ggplot_theme(methylation_change_boxplots, 
  ylab = "Mean Methylation Change", xlab = "Genomic Feature", show_legend = F, 
  scale_y = scale_y_continuous(limits = c(-0.5, 0.5)), 
  fill_colors = feature_colors[c(1:9, 12)])
methylation_change_boxplots = methylation_change_boxplots + 
  scale_x_discrete(labels = gsub("TMRs", "Clusters", levels(genomic_feature_mean_methylation_change_selected$region_type)))
ggsave(plot = methylation_change_boxplots, "../../presentation_figures/methylation_change_boxplots.pdf", width = 16, height = 9)
