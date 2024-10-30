# Create a function which will plot the CpG correlations for a specified transcript
plot_cpg_cor_values_transcript = function(transcript, title = NULL, 
  ylabel = "DNA Methylation-Transcription Correlation", xlabel = "Distance of CpG to TSS (bp)"){
  
  plot = methodical::plotMethSiteCorCoefs(meth_site_cor_values = mcrpc_correlation_results[[transcript]],
    reference_tss = plyranges::filter(tss_gr, transcript_id == transcript),
    xlabel = xlabel,  ylabel = ylabel, title = title) +
    geom_hline(yintercept = 0, linetype = "dashed") +
    scale_x_continuous(limits = c(-5000, 5000), breaks = seq(-4000, 4000, 2000), exp= c(0.005, 0.005), labels = scales::comma) +
    theme(axis.title.x = ggtext::element_markdown(hjust = 0.5, size = 20))
  
  return(plot)
  
}

# Create a function which will plot the promoter correlations for a specified transcript
plot_promoter_correlations = function(transcript_id, title = "Promoter\nMethylation Change"){
  
  # Get differential methylation results for transcript
  transcript_results = dplyr::filter(mcrpc_sample_correlation_tables_combined, table1_feature == transcript_id)
  
  # Create plot
  ggplot(transcript_results, 
    aes(y = cor, x = definition, fill = definition, label = significance)) + 
    geom_col(color = "black") + 
    geom_text(vjust = ifelse(transcript_results$cor >= 0, "bottom", "top"), size = 8) + 
    geom_hline(yintercept = 0) +
    scale_fill_manual(values = RColorBrewer::brewer.pal(9, "YlGn")[c(2, 4, 6, 8, 9)]) +
    theme_classic() +
    scale_y_continuous(exp= expansion(mult = c(0.05, 0.05))) +
    labs(x = "Promoter Definition", y = title, 
      title = NULL) +
    theme(plot.title = ggtext::element_markdown(hjust = 0.5, size = 18), 
      axis.title.y = ggtext::element_markdown(size = 20), axis.title.x = ggtext::element_markdown(size = 20),
      axis.text = element_text(size = 18), axis.text.x = element_text(hjust = 1), legend.position = "None")

}