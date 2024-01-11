# Find differentially methylated clusters in CPGEA samples for different cluster sources

# Load requried packages
library(dplyr)
library(ggplot2)
source("~/r_functions/table_tests.R")

# Get different sources of clusters
cpgea_normal_clusters = readRDS("../final_tmr_granges/cpgea_normal_tmrs_5kb.rds")
cpgea_tumour_clusters = readRDS("../final_tmr_granges/cpgea_tumour_tmrs_5kb.rds")
mcrpc_clusters = readRDS("../final_tmr_granges/mcrpc_tmrs_5kb.rds")

cpgea_normal_clusters$tmr_name = paste0("cpgea_normal_", cpgea_normal_clusters$tmr_name)
cpgea_tumour_clusters$tmr_name = paste0("cpgea_tumour_", cpgea_tumour_clusters$tmr_name)
mcrpc_clusters$tmr_name = paste0("mcrpc_", mcrpc_clusters$tmr_name)

# Combine TMRs
all_tmrs = c(cpgea_normal_clusters, cpgea_tumour_clusters, mcrpc_clusters)

# Load methylation tables
cpgea_normal_cluster_methylation = data.frame(data.table::fread("tmr_methylation_tables/cpgea_normal_tmrs_meth_all_samples.tsv.gz"), row.names = 1)
cpgea_tumour_cluster_methylation = data.frame(data.table::fread("tmr_methylation_tables/cpgea_tumour_tmrs_meth_all_samples.tsv.gz"), row.names = 1)
mcrpc_cluster_methylation = data.frame(data.table::fread("tmr_methylation_tables/mcrpc_tmrs_meth_all_samples.tsv.gz"), row.names = 1)

row.names(cpgea_normal_cluster_methylation) = paste0("cpgea_normal_", row.names(cpgea_normal_cluster_methylation))
row.names(cpgea_tumour_cluster_methylation) = paste0("cpgea_tumour_", row.names(cpgea_tumour_cluster_methylation))
row.names(mcrpc_cluster_methylation) = paste0("mcrpc_", row.names(mcrpc_cluster_methylation))

# Combine methylation tables
all_cluster_methylation_table = bind_rows(
  cpgea_normal = cpgea_normal_cluster_methylation,
  cpgea_tumour = cpgea_tumour_cluster_methylation,
  mcrpc = mcrpc_cluster_methylation, .id = "cluster_source"
)

all_cluster_methylation_table = all_cluster_methylation_table[all_tmrs$tmr_name, ]

all_cluster_methylation_change = select(all_cluster_methylation_table, starts_with("T")) - select(all_cluster_methylation_table, starts_with("N"))

# Calculate differential methylation using t-tests
system.time({all_cluster_methylation_change_results = 
  table_tests(test_table = all_cluster_methylation_change, test = t.test, gene_names = all_tmrs$gene_name)})

# Add cluster direction and methylation change directionto results
all_cluster_methylation_change_results$cluster_direction = all_tmrs$direction
all_cluster_methylation_change_results$meth_change_direction = ifelse(all_cluster_methylation_change_results$q_value > 0.05, "Unchanged", 
  ifelse(sign(all_cluster_methylation_change_results$mean_change) == 1, "Hypermethylated", "Hypomethylated"))
all_cluster_methylation_change_results$meth_change_direction = factor(all_cluster_methylation_change_results$meth_change_direction, levels = c("Hypomethylated", "Unchanged", "Hypermethylated"))

# Remove results with NA values
all_cluster_methylation_change_results = filter(all_cluster_methylation_change_results, !is.na(p_value))

#
all_cluster_methylation_change_results$cluster_source = gsub("_ENST.*", "", all_cluster_methylation_change_results$feature)

# Get number of hypermethylated and hypomethylated clusters for negative and positive clusters
all_cluster_direction_meth_change = mutate(summarize(group_by(all_cluster_methylation_change_results, 
  cluster_source, cluster_direction, meth_change_direction), count = dplyr::n()), proportion = count/sum(count)) 

# Make barplot for differential methylation of clusters for different cluster types
cluster_direction_meth_change_barplot = ggplot(all_cluster_direction_meth_change, 
  aes(x = cluster_source, y = count, fill = meth_change_direction, label = paste0(round(proportion, 2)*100, "%"))) +
  geom_col(position = "dodge", color = "black") + 
  geom_text(mapping = aes(x = cluster_source, y = count + 100, group = meth_change_direction), position = position_dodge(width = 0.9), size = 4) +
  theme_classic() +
  theme(plot.title = element_text(hjust = 0.5, size = 24), 
  	axis.title = element_text(size = 20), axis.text = element_text(size = 18), 
  	legend.text = element_text(size = 18), legend.title = element_text(size = 20), legend.key.size = unit(1, "cm"),
    strip.text.x = element_text(size = 20)) +
  scale_x_discrete(labels = c("Normal\nProstate", "Prostate\nTumours", "Prostate\nMetastases")) + 
  scale_y_continuous(expand = expansion(mult = c(0, 0.05))) + 
  scale_fill_manual(values = c("#4B878BFF", "grey", "#D01C1FFF")) +
  labs(x = "TMR Source", y = "Number of TMRs", fill = NULL,  title = NULL) +
  facet_wrap(as.formula(paste("~", "cluster_direction")), nrow = 2, 
    labeller = as_labeller(setNames(c("Negative TMRs", "Positive TMRs"), c("Negative", "Positive"))))
cluster_direction_meth_change_barplot = ggpubr::ggarrange(cluster_direction_meth_change_barplot, labels = "A")
cluster_direction_meth_change_barplot
ggsave(plot = cluster_direction_meth_change_barplot, "tmr_differential_meth_barplot.pdf", width = 16, height = 9)
ggsave(plot = cluster_direction_meth_change_barplot, "tmr_differential_meth_barplot_1.pdf", width = 8, height = 9)

# Get the hypermethylated negative TMRs for each dataset
hypermethylated_negative_tmrs = with(filter(all_cluster_methylation_change_results, 
  cluster_direction == "Negative" & meth_change_direction == "Hypermethylated"), 
  lapply(split(gene, cluster_source), unique))

# Get the hypomethylated negative TMRs for each dataset
hypomethylated_negative_tmrs = with(filter(all_cluster_methylation_change_results, 
  cluster_direction == "Negative" & meth_change_direction == "Hypomethylated"), 
  lapply(split(gene, cluster_source), unique))

# Get the hypermethylated positive TMRs for each dataset
hypermethylated_positive_tmrs = with(filter(all_cluster_methylation_change_results, 
  cluster_direction == "Positive" & meth_change_direction == "Hypermethylated"), 
  lapply(split(gene, cluster_source), unique))

# Get the hypomethylated positive TMRs for each dataset
hypomethylated_positive_tmrs = with(filter(all_cluster_methylation_change_results, 
  cluster_direction == "Positive" & meth_change_direction == "Hypomethylated"), 
  lapply(split(gene, cluster_source), unique))

# Save results as a list
tmr_diff_meth_genes = list(
  hypermethylated_negative = hypermethylated_negative_tmrs,
  hypomethylated_negative = hypomethylated_negative_tmrs,
  hypermethylated_positive = hypermethylated_positive_tmrs,
  hypomethylated_positive = hypomethylated_positive_tmrs
  )
saveRDS(tmr_diff_meth_genes, "tmr_diff_meth_genes.rds")


