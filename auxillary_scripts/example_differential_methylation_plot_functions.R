# Create a function which plots promoter methylation change using the different promoter definitions
plot_promoter_methylation_change = function(transcript){
  
  # Get gene name from transcript
  gene = gencode_tss_gr$gene_name[gencode_tss_gr$transcript_id == transcript]
  
  # Get differential methylation results for transcript
  transcript_results = dplyr::filter(promoter_diff_meth_results, transcript_id == transcript)
  
  # Create plot
  ggplot(transcript_results, 
    aes(y = meth_diff, x = definition, fill = definition, label = significance)) + 
    geom_col(color = "black") + 
    geom_text(vjust = ifelse(transcript_results$meth_diff >= 0, "bottom", "top"), size = 8) + 
    geom_hline(yintercept = 0) +
    scale_fill_manual(values = RColorBrewer::brewer.pal(9, "YlGn")[c(2, 4, 6, 8, 9)]) +
    theme_classic() +
    scale_y_continuous(expand = expansion(mult = c(0.05, 0.05))) +
    labs(x = "Promoter Definition", y = paste(paste0("*", gene, "*"), "Promoter<br>Methylation Change"), 
      title = NULL) +
    theme(plot.title = ggtext::element_markdown(hjust = 0.5, size = 20), 
      axis.title.y = ggtext::element_markdown(size = 20), axis.title.x = element_text(size = 20),
      axis.text = element_text(size = 18), axis.text.x = element_text(hjust = 1), legend.position = "None") 

}

# Make a function which will plot CpG methylation change for all CpGs in +/- 5KB of a TSS
plot_cpg_methylation_change = function(transcript, title = NULL){
  
  # Get gene name from transcript
  gene = gencode_tss_gr$gene_name[gencode_tss_gr$transcript_id == transcript]
  
  # Get TSS of transcript
  tss_gr = plyranges::filter(gencode_tss_gr, transcript_id == transcript)
  
  # Get methylation values of all CpG values within +/- 5 KB of the TSS
  transcript_methylation = methodical::extractGRangesMethSiteValues(cpgea_meth_rse, promoters(tss_gr, 5000, 5001))
  
  # Get mean methylation change for each site
  transcript_methylation_change = data.frame(meth_change = 
      rowMeans(select(transcript_methylation, starts_with("T")) 
        - select(transcript_methylation, starts_with("N")), na.rm = T))
  
  # Plot methylation change for the CpGs
  cpg_plot = methodical::plotMethylationValues(meth_site_values = transcript_methylation_change, 
    sample_name = "meth_change", reference_tss = tss_gr, 
    xlabel = paste("Distance to", paste0("*", gene, "*"), "TSS (bp)"), ylabel = "DNA Methylation Change", title = title) +
    geom_hline(yintercept = 0, linetype = "dashed") +
    scale_x_continuous(limits = c(-5000, 5000), breaks = seq(-4000, 4000, 2000), expand = c(0.005, 0.005), labels = scales::comma) +
    scale_y_continuous(breaks = seq(-1, 1, 0.2), expand = c(0.05, 0.05), labels = scales::comma) +
    theme(axis.title.x = ggtext::element_markdown(size = 18))
  
}