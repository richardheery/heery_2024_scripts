# Evaluate TMRs associated with MANE transcripts in TCGA 

# Load required packages
library(methodical)
library(doParallel)
source("../auxillary_scripts/plotting_functions.R")

# Load the names of MANE transcripts
mane_transcripts = readRDS("../auxillary_data/genomic_annotation/mane_transcript_ids.rds")

# Get CPGEA tumour TMRs and filter for those associated with MANE transcripts
mcrpc_tmrs = readRDS("../finding_tmrs/mcrpc_tmrs.rds")
mcrpc_tmrs = mcrpc_tmrs[mcrpc_tmrs$ID %in% mane_transcripts]
mcrpc_tmrs$gene_id = gsub("\\.[0-9]", "", mcrpc_tmrs$gene_id)

# Rename TMRs with gene ID
mcrpc_tmrs$tmr_name = paste0(mcrpc_tmrs$gene_id, gsub(".*_tmr_", "_tmr_", mcrpc_tmrs$tmr_name))

# Get paths to DESeq2 normalized count tables for TCGA projects
deseq2_normalized_count_files = list.files("../auxillary_data/rnaseq_data/tcga_deseq_normalized_counts", full.names = T)
names(deseq2_normalized_count_files) = gsub("_deseq_normalized_counts.tsv.gz", "", basename(deseq2_normalized_count_files))

# Load meth RSE for TCGA hg19
tcga_meth_rse_hg19 = HDF5Array::loadHDF5SummarizedExperiment("../auxillary_data/methylation_data/tcga_meth_array_hg19")
hg19tohg38_chain = rtracklayer::import.chain("../auxillary_data/genomic_annotation/hg19ToHg38.over.chain")
hg38_cpgs = extractMethSitesFromGenome("BSgenome.Hsapiens.UCSC.hg38")
tcga_meth_rse_hg38 = liftoverMethRSE(tcga_meth_rse_hg19, chain = hg19tohg38_chain, permitted_target_regions = hg38_cpgs)

cl = parallel::makeCluster(11)
doParallel::registerDoParallel(cl, 11)

# Calculate correlation for CPGEA tumour TMRs in tumour samples for each project
# Took 22 minutes with 3 cores
system.time({mcrpc_tmr_tcga_tumour_cors = foreach(project = names(deseq2_normalized_count_files))  %dopar% { 
  
  # Print name of current project
  message(paste("Starting", project))
  
  # Read in counts table for project
  counts_table = data.frame(data.table::fread(deseq2_normalized_count_files[project]), row.names = 1)
  
  # Find common samples to tcga_meth_rse_hg38 and counts table
  common_samples = intersect(names(counts_table), colnames(tcga_meth_rse_hg38))
  
  # Subset for tumour samples
  common_tumour_samples = common_samples[endsWith(common_samples, "01")]
  
  # If there are less than 50 tumour samples, skip to next iteration
  if(length(common_tumour_samples) < 50){
    message("Less than 50 common tumour samples so skipping to next project")
    return(NULL)
  }
  
  # Calculate correlation results for tumour samples
  cor_results = methodical::calculateRegionMethylationTranscriptCors(
    meth_rse = tcga_meth_rse_hg38, 
    transcript_expression_table = counts_table, genomic_regions = mcrpc_tmrs, 
    genomic_region_names = mcrpc_tmrs$tmr_name, 
    genomic_region_transcripts = mcrpc_tmrs$gene_id, samples_subset = common_tumour_samples, 
    cor_method = "s", BPPARAM = BiocParallel::MulticoreParam(workers = 20))
  
  # Return results
  cor_results
  
}})
names(mcrpc_tmr_tcga_tumour_cors) = names(deseq2_normalized_count_files)
saveRDS(mcrpc_tmr_tcga_tumour_cors, "mcrpc_tmr_tcga_tumour_cors.rds")
mcrpc_tmr_tcga_tumour_cors = readRDS("mcrpc_tmr_tcga_tumour_cors.rds")

# Combine results for different cancer types
mcrpc_tmr_tcga_tumour_cors = dplyr::bind_rows(mcrpc_tmr_tcga_tumour_cors, .id = "project")

# Add direction of TMRs 
mcrpc_tmr_tcga_tumour_cors$direction = mcrpc_tmrs$direction[match(mcrpc_tmr_tcga_tumour_cors$genomic_region_name, mcrpc_tmrs$tmr_name)]

# Divide the cancer types into 3 groups
mcrpc_tmr_tcga_tumour_cors$group = ceiling(as.numeric(as.factor(mcrpc_tmr_tcga_tumour_cors$project))/10)
mcrpc_tmr_tcga_tumour_cors$project = factor(mcrpc_tmr_tcga_tumour_cors$project)
mcrpc_tmr_tcga_tumour_cors$cancer_type = mcrpc_tmr_tcga_tumour_cors$project
levels(mcrpc_tmr_tcga_tumour_cors$cancer_type) = c("Adrenocortical", "Bladder", "Breast", "Cervical", "Colon", "Esophageal", "Glioma", "Head and\nNeck", 
    "Kidney\nChromophobe", "Kidney\nClear Cell", "Kidney\nPapillary Cell", "Lower Grade\nGlioma", "Liver", "Lung Adeno-\ncarcinoma", "Lung\nSquamous", "Mesothelioma",
    "Pancreatic", "PCPG", "Prostate", "Rectal", "Sarcoma", "Skin\nMelanoma", "Stomach",
    "Testicular\nGerm Cell", "Thyroid", "Thymoma", "Uterine\nCarcinosarcoma", "Endometrial", "Uveal\nMelanoma")

# Create boxplots of the correlations
mcrpc_tmr_tcga_tumour_cors_boxplots = ggplot(mcrpc_tmr_tcga_tumour_cors, aes(x = cancer_type, y = cor, fill = direction)) +
  geom_boxplot() + geom_hline(yintercept = 0, linetype = "dashed") +
  facet_wrap("~ group", nrow = 3, scales = "free")
mcrpc_tmr_tcga_tumour_cors_boxplots = customize_ggplot_theme(mcrpc_tmr_tcga_tumour_cors_boxplots, fill_colors = colour_list$purple_and_gold_light, 
  xlab = "Cancer Type", ylab = "TMR Methylation-Expression Correlation", fill_title = "Direction", axis_text_size = 16) +
  theme(strip.background = element_blank(), strip.text = element_blank())
mcrpc_tmr_tcga_tumour_cors_boxplots = ggpubr::ggarrange(mcrpc_tmr_tcga_tumour_cors_boxplots, labels = "D")
mcrpc_tmr_tcga_tumour_cors_boxplots
ggsave(plot = mcrpc_tmr_tcga_tumour_cors_boxplots, "../figures/figure6D.pdf", width = 16, height = 9)
