# Calculate promoter methylation in CPGEA samples for different promoter definitions 

# Load required packages
library(methodical)
library(dplyr)
library(plotR)

# Get promoter definition list
promoter_definition_list = readRDS("promoter_definition_list.rds")
promoter_definition_list = lapply(promoter_definition_list, function(x) setNames(x, x$transcript_id))

# Get path to CPGEA methrix and convert a methylation RSE
cpgea_rse = HDF5Array::loadHDF5SummarizedExperiment("../auxiallary_data/methylation_data/cpgea_meth_rse_hg38")

# Create a BPPARAM object
bpparam = BiocParallel::MulticoreParam(3)

# Make a directory to store result files
mkdir("promoter_definition_methylation_tables")

# Get methylation values for all promoter methylation definitions and save to promoter_definition_methylation_tables. Took 20 minutes. 
system.time({for(definition in names(promoter_definition_list)){
  message(definition)
  promoter_definition_methylation = methodical::summarizeRegionMethylation(
    meth_rse = cpgea_rse, genomic_regions = promoter_definition_list[[definition]], 
    BPPARAM = bpparam)
  promoter_definition_methylation = tibble::rownames_to_column(promoter_definition_methylation, "transcript_id")
  data.table::fwrite(promoter_definition_methylation, 
    file = sprintf("promoter_definition_methylation_tables/%s_definition_methylation_table.tsv.gz", definition), 
    sep = "\t", row.names = F, na = "NA", quote = F)
}})

# Get paths to all promoter methylation definition tables
promoter_definition_methylation_tables = list.files("promoter_definition_methylation_tables", full.names = T, pattern = "_methylation_table.tsv.gz")
names(promoter_definition_methylation_tables) = gsub("_definition_methylation_table.tsv.gz", "", basename(promoter_definition_methylation_tables))

# Read in tables as a list and save
promoter_definition_methylation_tables = lapply(promoter_definition_methylation_tables, function(x)
  data.frame(data.table::fread(x)))
for(definitiom in names(promoter_definition_methylation_tables)){
  row.names(promoter_definition_methylation_tables[[definition]]) = cage_transcripts
}

# Save list of promoter methylation tables
saveRDS(promoter_definition_methylation_tables, "promoter_definition_methylation_tables/promoter_definition_methylation_tables.rds")
promoter_definition_methylation_tables = readRDS("promoter_definition_methylation_tables/promoter_definition_methylation_tables.rds")