# Download methylation and RNA-seq data from TumourMethData

# Load required packages
library(TumourMethData)

# Download WGBS data from CPGEA, MCRPC and TCGA along with TCGA array data
download_meth_dataset(dataset = "cpgea_wgbs_hg38", "methylation_data/")
download_meth_dataset(dataset = "mcrpc_wgbs_hg38", "methylation_data/")
download_meth_dataset(dataset = "tcga_wgbs_hg38", "methylation_data/")
download_meth_dataset(dataset = "tcga_450k_array_hg19", "methylation_data/")
download_meth_dataset(dataset = "roadmap_wgbs_hg38", "methylation_data/")


# Download transcript counts for CPGEA and MCRPC 
download_rnaseq_dataset(dataset = "cpgea_wgbs_hg38", "rnaseq_data/")
download_rnaseq_dataset(dataset = "mcrpc_wgbs_hg38", "rnaseq_data/")